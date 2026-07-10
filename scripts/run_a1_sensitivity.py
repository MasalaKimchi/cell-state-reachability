"""A1 — verdict sensitivity radius, calibrated in units of measured standard error.

The reachability verdict is a point estimate on noisy effect vectors. This asks:
how far would the effects have to be wrong to flip it? We answer with two arms:

  (M) MEASUREMENT-ERROR robustness — the effect estimates carry a measured standard
      error (atlas_lfcSE.npz). Resample E ~ N(E_hat, SE^2), re-solve, and report the
      bootstrap distribution of the reachable cosine. Because undirected noise ADDS
      anisotropy, it tends to *inflate* the reachable cosine (reachability-by-chance),
      so the bootstrap LOWER bound vs the null threshold is the honest robustness number.

  (C) COORDINATED-BIAS radius — the worst-case adversary attenuates the target-aligned
      component of every generator by a fraction beta (E' = E - beta * (E d_hat) d_hat^T
      on the readout). This is re-solve-proof (it removes the d-direction from the
      dictionary) and is the optimal attack; we bisect for the beta* that drops the
      cosine to the null p99, and express the induced per-entry shift in units of the
      median measured SE. A near-threshold (graded-reachable) verdict has a small beta*,
      which is *why* the design-based randomization argument is load-bearing.

Runs one culture condition (default Stim48hr, the differentiated state) across all four
target axes. Memory-frugal: loads E and SE for the condition once (float32), reuses a
single working buffer, never materializes a float64 full copy.
"""
import numpy as np, json, os, sys, time, gc
import os as _os
_REPO_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
_os.chdir(_REPO_ROOT); sys.path.insert(0, _REPO_ROOT)
import reachability as R

COND = sys.argv[1] if len(sys.argv) > 1 else "Stim48hr"
N_BOOT = int(os.environ.get("A1_NBOOT", "12"))
N_PROBE = int(os.environ.get("A1_NPROBE", "8"))
TARGETS = ["toward_Th1", "toward_Th2", "toward_younger", "toward_older"]
WORK = "analysis_cache/atlas_work"
ART_INPUTS = os.environ["A1_INPUTS"]   # path to inputs.npz
ART_SE = os.environ["A1_SE"]           # path to atlas_lfcSE.npz

ip = np.load(ART_INPUTS, allow_pickle=True, mmap_mode="r")
se = np.load(ART_SE, allow_pickle=True, mmap_mode="r")
E_c = np.asarray(ip[f"E_{COND}"])          # (G, V) float32  ~283MB
SE_c = np.asarray(se[f"SE_{COND}"])        # (G, V) float32  ~283MB
assert E_c.shape == SE_c.shape
gene_c = ip[f"gene_{COND}"]
print(f"loaded {COND}: E {E_c.shape}, SE {SE_c.shape}", flush=True)

_buf = np.empty_like(E_c)                  # single reused working buffer

def solve_cos(mat, d, mask):
    return float(R.reachability(mat, d, hvg_mask=mask).reachable_cosine)

def boot_draw(seed, d, mask):
    rng = np.random.default_rng(seed)
    rng.standard_normal(size=E_c.shape, dtype=np.float32, out=_buf)
    np.multiply(_buf, SE_c, out=_buf)      # k=1: draw at the measured SE
    np.add(_buf, E_c, out=_buf)            # E_hat + SE*z
    return solve_cos(_buf, d, mask)

def cos_attenuated(beta, d, mask, dhat, proj):
    np.copyto(_buf, E_c)
    # subtract beta * proj_g * dhat over masked columns
    _buf[:, mask] -= (beta * np.outer(proj, dhat)).astype(np.float32)
    return solve_cos(_buf, d, mask)

results = []
for name in TARGETS:
    t0 = time.time()
    d = np.asarray(ip[f"t_{name}"]); mask = d != 0
    base = solve_cos(E_c, d, mask)
    # null threshold (analytic anisotropy null): p99 = mean + 2.326 sd
    nl = R.analytic_anisotropy_null(E_c, d, hvg_mask=mask, n_probe=N_PROBE, observed=base)
    p99 = float(nl.null_mean + 2.326 * nl.null_sd)
    z = float(nl.z)
    # (M) measurement-error bootstrap
    boot = [boot_draw(1000 + i, d, mask) for i in range(N_BOOT)]
    boot = np.array(boot)
    lo, hi = float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))
    # (C) coordinated-bias radius: bisection on beta in [0,1] to hit p99
    dm = d[mask].astype(np.float64); dhat = dm / np.linalg.norm(dm)
    proj = (E_c[:, mask].astype(np.float64) @ dhat)   # (G,) transient float64 vector, not matrix
    med_se = float(np.median(SE_c[:, mask]))
    def kse(beta):
        shift = beta * np.abs(np.outer(proj, dhat))
        v = float(np.sqrt((shift**2).mean()) / med_se); del shift; return v
    blo, bhi = 0.0, 1.0
    c_hi = cos_attenuated(bhi, d, mask, dhat, proj)
    if c_hi > p99:
        beta_star, kse_star = float("nan"), float("nan")   # even full removal doesn't flip
    else:
        for _ in range(8):
            bm = 0.5 * (blo + bhi)
            if cos_attenuated(bm, d, mask, dhat, proj) > p99:
                blo = bm
            else:
                bhi = bm
        beta_star = 0.5 * (blo + bhi); kse_star = kse(beta_star)
    row = dict(target=name.replace("toward_", ""), condition=COND,
               baseline_cosine=base, null_p99=p99, null_z=z,
               boot_lo=lo, boot_hi=hi, boot_mean=float(boot.mean()),
               boot_clears_null=bool(lo > p99), n_boot=N_BOOT,
               beta_star=beta_star, kse_star=kse_star, med_se=med_se,
               seconds=round(time.time() - t0, 1))
    results.append(row)
    print(f"[{name}] base={base:.4f} p99={p99:.4f} z={z:.2f} "
          f"boot95=[{lo:.4f},{hi:.4f}] clears_null={row['boot_clears_null']} "
          f"beta*={beta_star:.3f} kSE*={kse_star:.3f} ({row['seconds']}s)", flush=True)
    gc.collect()

import pandas as pd
df = pd.DataFrame(results)
df.to_csv(f"{WORK}/a1_sensitivity_radius.csv", index=False)
df.to_csv("results/a1_sensitivity_radius.csv", index=False)
print("A1 COMPLETE\n" + df.to_string(index=False), flush=True)
