"""Extract per-donor single-gene effect vectors for leakage-free LODO.

Streams only the planned pseudobulk rows from the 44.6 GB remote AnnData
(GWCD4i.pseudobulk_merged.h5ad) over HTTP range requests, normalizes each
pseudobulk to log1p-CPM, and builds, PER DONOR, a single-gene effect matrix
    atoms_d[g] = mean_over_guides( log1pCPM(guide g, donor d) ) - mean( log1pCPM(NTC, donor d) )
No cross-donor pooling: donor d's atoms use ONLY donor d's cells. That is what
makes leave-one-donor-out leakage-free.

Output: repo/results/tcell_lodo_substrate.npz
Progress: repo/_extract_progress.json (updated each row).
"""
import io, json, time, sys, os
import numpy as np
import pandas as pd
import requests
import h5py

HERE = os.path.dirname(os.path.abspath(__file__))
URL = ("https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/"
       "marson2025_data/GWCD4i.pseudobulk_merged.h5ad")
PROG = os.path.join(HERE, "_extract_progress.json")
OUT  = os.path.join(HERE, "results", "tcell_lodo_substrate.npz")


class HTTPRangeFile(io.RawIOBase):
    """Seekable read-only file over HTTP range requests with a small block cache."""
    def __init__(self, url, block=1024 * 1024, session=None):
        self.url, self.block = url, block
        self.sess = session or requests.Session()
        self._pos = 0
        self._cache = {}
        self.n_requests = 0
        self.bytes_fetched = 0
        r = self.sess.head(url, timeout=30); r.raise_for_status()
        self.size = int(r.headers["Content-Length"])

    def _fetch_block(self, bi):
        if bi in self._cache:
            return self._cache[bi]
        start = bi * self.block
        end = min(start + self.block, self.size) - 1
        r = self.sess.get(self.url, headers={"Range": f"bytes={start}-{end}"}, timeout=120)
        r.raise_for_status()
        data = r.content
        self.n_requests += 1
        self.bytes_fetched += len(data)
        if len(self._cache) > 80:
            self._cache.pop(next(iter(self._cache)))
        self._cache[bi] = data
        return data

    def readable(self): return True
    def seekable(self): return True
    def seek(self, off, whence=io.SEEK_SET):
        if whence == io.SEEK_SET: self._pos = off
        elif whence == io.SEEK_CUR: self._pos += off
        elif whence == io.SEEK_END: self._pos = self.size + off
        return self._pos
    def tell(self): return self._pos
    def read(self, n=-1):
        if n is None or n < 0: n = self.size - self._pos
        n = min(n, self.size - self._pos)
        if n <= 0: return b""
        out = bytearray(); pos = self._pos
        while len(out) < n:
            bi = pos // self.block
            off = pos - bi * self.block
            blk = self._fetch_block(bi)
            take = min(len(blk) - off, n - len(out))
            out += blk[off:off + take]; pos += take
        self._pos = pos
        return bytes(out)
    def readinto(self, b):
        d = self.read(len(b)); b[:len(d)] = d; return len(d)


def main():
    t0 = time.time()
    atom_plan = pd.read_csv(os.path.join(HERE, "_lodo_atom_plan.csv"))
    ntc_plan = pd.read_csv(os.path.join(HERE, "_lodo_ntc_plan.csv"))
    axis = np.load(os.path.join(HERE, "_lodo_axis.npz"), allow_pickle=True)
    var_gene_name = axis["var_gene_name"]
    panel_genes = list(axis["panel_genes"])
    NG = len(var_gene_name)
    donors = ["D1", "D2", "D3", "D4"]

    rf = HTTPRangeFile(URL, block=1024 * 1024)
    f = h5py.File(rf, "r", driver="fileobj",
                  rdcc_nbytes=256 * 1024 * 1024, rdcc_nslots=100003)
    indptr = f["X"]["indptr"][:]
    data_ds = f["X"]["data"]; idx_ds = f["X"]["indices"]

    def read_logcpm(row):
        lo, hi = int(indptr[row]), int(indptr[row + 1])
        cols = idx_ds[lo:hi]; vals = data_ds[lo:hi]
        v = np.zeros(NG, dtype=np.float64)
        v[cols] = vals
        tot = v.sum()
        if tot <= 0:
            return v
        return np.log1p(v / tot * 1e4)  # log1p-CP10K

    all_rows = list(atom_plan["row"]) + list(ntc_plan["row"])
    n_total = len(all_rows)
    logcpm = {}
    done = 0
    for row in all_rows:
        logcpm[int(row)] = read_logcpm(int(row))
        done += 1
        if done % 5 == 0 or done == n_total:
            json.dump({"done": done, "total": n_total,
                       "elapsed_s": round(time.time() - t0, 1),
                       "http_req": rf.n_requests,
                       "MB": round(rf.bytes_fetched / 1e6, 1)},
                      open(PROG, "w"))

    # Per-donor NTC baseline (mean log-CPM over that donor's NTC pseudobulks)
    ntc_baseline = {}
    for d in donors:
        rows = ntc_plan[ntc_plan.donor == d]["row"].astype(int).tolist()
        ntc_baseline[d] = np.mean([logcpm[r] for r in rows], axis=0)

    # Per-donor single-gene atoms: mean over that gene's guides in that donor, minus donor NTC
    atoms = {d: np.zeros((len(panel_genes), NG)) for d in donors}
    coverage = {d: np.zeros(len(panel_genes), dtype=bool) for d in donors}
    for gi, g in enumerate(panel_genes):
        for d in donors:
            rows = atom_plan[(atom_plan.perturbed_gene_name == g) &
                             (atom_plan.donor == d)]["row"].astype(int).tolist()
            if not rows:
                continue
            mean_g = np.mean([logcpm[r] for r in rows], axis=0)
            atoms[d][gi] = mean_g - ntc_baseline[d]
            coverage[d][gi] = True

    np.savez_compressed(
        OUT,
        panel_genes=np.array(panel_genes),
        var_gene_name=var_gene_name,
        atoms_D1=atoms["D1"], atoms_D2=atoms["D2"],
        atoms_D3=atoms["D3"], atoms_D4=atoms["D4"],
        coverage_D1=coverage["D1"], coverage_D2=coverage["D2"],
        coverage_D3=coverage["D3"], coverage_D4=coverage["D4"],
        ntc_D1=ntc_baseline["D1"], ntc_D2=ntc_baseline["D2"],
        ntc_D3=ntc_baseline["D3"], ntc_D4=ntc_baseline["D4"],
    )
    json.dump({"done": n_total, "total": n_total, "status": "COMPLETE",
               "elapsed_s": round(time.time() - t0, 1),
               "http_req": rf.n_requests, "MB": round(rf.bytes_fetched / 1e6, 1),
               "out": OUT}, open(PROG, "w"))
    print("DONE", OUT, "elapsed", round(time.time() - t0, 1), "s")


if __name__ == "__main__":
    main()
