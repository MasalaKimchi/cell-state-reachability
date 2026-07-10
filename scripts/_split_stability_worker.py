"""Certificate gene-axis split-stability test (held-out-modality analogue).

The Th2->Th1 target d splits into a knockdown-reachable part and an UP part
(the certificate's activation demand = positive residual after the NNLS cone
fit). We split the GENE axis into two random halves, fit the non-negative cone
on each half INDEPENDENTLY, and read the certificate (unmet upward demand) off
the full target axis from each half-fit. Stability of the certificate ranking
across the independent halves is the in-silico analogue of a held-out-modality
test: it asks whether the SAME activation genes surface no matter which half of
the transcriptome was used to build the knockdown cone.

Runs S seeds; each seed = two independent half-fits. Parallel across seeds.
"""
import numpy as np, json, time, sys, os
from scipy.optimize import nnls
from concurrent.futures import ProcessPoolExecutor

REPO = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "analysis_cache")
Zt = np.load(os.path.join(REPO, "split_inputs.npz"), allow_pickle=True)
E_m = Zt["E_m"].astype(np.float64)      # (P, G) generators on target-support genes
d_m = Zt["d_m"].astype(np.float64)      # (G,) target on those genes
names_m = Zt["names_m"]
G = d_m.shape[0]

def cert_score(rho, d):
    """unmet UPWARD demand: residual>0 AND target wants gene up."""
    return np.where((rho > 0) & (d > 0), rho, 0.0)

def one_seed(seed):
    rng = np.random.default_rng(seed)
    perm = rng.permutation(G)
    h1, h2 = perm[:G // 2], perm[G // 2:]
    t0 = time.time()
    w1, _ = nnls(E_m[:, h1].T, d_m[h1])     # fit cone on half 1 only
    w2, _ = nnls(E_m[:, h2].T, d_m[h2])     # fit cone on half 2 only
    rho1 = d_m - E_m.T @ w1                  # residual read on FULL target axis
    rho2 = d_m - E_m.T @ w2
    return dict(seed=int(seed),
                s1=cert_score(rho1, d_m).astype(np.float32),
                s2=cert_score(rho2, d_m).astype(np.float32),
                sec=round(time.time() - t0, 1))

if __name__ == "__main__":
    S = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    seeds = list(range(S))
    t0 = time.time()
    results = []
    for s in seeds:                        # serial (sandbox blocks ProcessPool semaphores)
        r = one_seed(s)
        results.append(r)
        print(f"seed {r['seed']} done ({r['sec']}s) [{len(results)}/{S}]", flush=True)
    S1 = np.vstack([r["s1"] for r in results])   # (S, G)
    S2 = np.vstack([r["s2"] for r in results])
    np.savez(os.path.join(REPO, "split_stability_raw.npz"),
             S1=S1, S2=S2, names_m=names_m, d_m=d_m, seeds=np.array(seeds))
    print(f"SAVED split_stability_raw.npz  S={S}  total {round(time.time()-t0,1)}s", flush=True)
