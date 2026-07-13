"""DEG-weighted evaluation of cell-state reachability — the Needles-in-the-Haystack
calibration applied to the reachability verdict.

Computes, for the headline targets:
  1. unweighted vs DEG-weighted reachable cosine + held-out-gene z (robustness recompute)
  2. positive-control ceiling (interpolated duplicate) + shuffled-target floor
  3. dynamic-range fraction = (observed - floor) / (ceiling - floor)

Datasets: Norman et al. 2019 (K562 CRISPRa, held-out perturbations) and the Tier-2 CD4
CRISPRi atlas (toward_Th1 across Rest/Stim8hr/Stim48hr). Writes tidy CSVs to results/.

Run:  python run_deg_weighted_eval.py [--quick]
  --quick  : Norman only + fewer shuffles (seconds); default also does Tier-2 (~25 min).
"""
import os, sys, json, time, argparse
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import reachability as R

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(ROOT, "notebooks", "cache")
DATA = os.path.join(ROOT, "data")
RESULTS = os.path.join(ROOT, "results")
SCHEME = "abs"          # DEG weight scheme: w_j = |d_j| (continuous effect-size weighting)
N_HELDOUT = 60          # shuffles for held-out-gene null (Norman; Tier-2 uses N_HELDOUT_T2)
N_FLOOR = 200           # shuffles for the negative-control (shuffled-target) floor (Norman)
N_CEIL = 60             # interpolated-duplicate positive-control targets (Norman)
# Tier-2 is ~18 s per full NNLS (6871 generators x 6188 genes), so Monte-Carlo counts are
# reduced there: the null MEAN and ceiling MEDIAN converge in a few dozen samples. Full
# floor+ceiling+DRF calibration is done on the HEADLINE Rest condition; the cheaper
# robustness recompute (cosine + held-out z, unweighted vs DEG-weighted) runs on all three.
N_HELDOUT_T2 = 40
N_FLOOR_T2 = 40
N_CEIL_T2 = 40
SEED = 0
# Gram-reuse + process-parallel fast path in reachability.py. -1 = all cores, 1 = serial.
# Accelerates the UNWEIGHTED nulls; weighted nulls take n_jobs as a harmless no-op (their
# per-iter matrix is rescaled by sqrt(weights) so it is not a fixed dictionary). Output is
# identical to serial for both. Override with REACH_N_JOBS=1 to force serial.
N_JOBS = int(os.environ.get("REACH_N_JOBS", "-1"))


def calib_row(E, d, *, dataset, target, weights_scheme=None, n_floor=N_FLOOR,
              n_ceil=N_CEIL, n_heldout=N_HELDOUT, seed=SEED, extra=None):
    """One fully-calibrated verdict row (unweighted OR weighted, per weights_scheme)."""
    wv = R.deg_weights(d, scheme=weights_scheme) if weights_scheme else None
    res = R.reachability(E, d, weights=wv)
    ho = R.held_out_gene_validation(E, d, weights=wv, n_shuffles=n_heldout, seed=seed, n_jobs=N_JOBS)
    floor = R.shuffled_target_null(E, d, n_iter=n_floor, seed=seed, weights=wv, n_jobs=N_JOBS)
    pc = R.positive_control_ceiling(E, weights_scheme=weights_scheme, n_targets=n_ceil,
                                    support_size=5, noise=0.25, seed=seed)
    cal = R.calibrate_reachability(res.reachable_cosine, floor.null_cosines.mean(),
                                   pc["ceiling"], weighted=wv is not None,
                                   n_ceiling=pc["n_targets"])
    row = dict(dataset=dataset, target=target,
               weighting=("DEG-weighted(|d|)" if weights_scheme else "unweighted"),
               n_generators=int(E.shape[0]), n_genes=int(E.shape[1]),
               reachable_cosine=round(float(res.reachable_cosine), 4),
               residual_norm=round(float(res.residual_norm), 4),
               heldout_cosine=round(float(ho.held_out_cosine), 4),
               heldout_z=round(float(ho.z), 2),
               floor_mean=round(float(floor.null_cosines.mean()), 4),
               floor_p95=round(float(floor.p95), 4),
               ceiling=round(float(pc["ceiling"]), 4),
               ceiling_p25=round(float(pc["ceiling_p25"]), 4),
               ceiling_p75=round(float(pc["ceiling_p75"]), 4),
               dynamic_range_fraction=round(float(cal.dynamic_range_fraction), 4),
               support_size=int(res.support.size))
    if extra:
        row.update(extra)
    return row


def robust_row(E, d, *, dataset, target, weights_scheme=None, n_heldout=N_HELDOUT,
               seed=SEED, extra=None):
    """Lighter row: reachable cosine + held-out-gene z only (no floor/ceiling/DRF).
    Used where the full Monte-Carlo calibration is too expensive to repeat per condition."""
    wv = R.deg_weights(d, scheme=weights_scheme) if weights_scheme else None
    res = R.reachability(E, d, weights=wv)
    ho = R.held_out_gene_validation(E, d, weights=wv, n_shuffles=n_heldout, seed=seed, n_jobs=N_JOBS)
    row = dict(dataset=dataset, target=target,
               weighting=("DEG-weighted(|d|)" if weights_scheme else "unweighted"),
               n_generators=int(E.shape[0]), n_genes=int(E.shape[1]),
               reachable_cosine=round(float(res.reachable_cosine), 4),
               residual_norm=round(float(res.residual_norm), 4),
               heldout_cosine=round(float(ho.held_out_cosine), 4),
               heldout_z=round(float(ho.z), 2),
               floor_mean=np.nan, floor_p95=np.nan, ceiling=np.nan,
               ceiling_p25=np.nan, ceiling_p75=np.nan, dynamic_range_fraction=np.nan,
               support_size=int(res.support.size))
    if extra:
        row.update(extra)
    return row


def load_norman():
    z = np.load(os.path.join(CACHE, "norman_effect_bundle.npz"), allow_pickle=True)
    return (z["Es"].astype(np.float64), z["Ed"].astype(np.float64),
            z["sing_labels"].astype(str), z["doub_labels"].astype(str),
            z["comp_idx"], z["comp_ok"].astype(bool))


def norman_targets(Es, Ed, sing, doub, comp_idx, comp_ok, k=6):
    """Headline Norman targets: held-out DOUBLE perturbations fit by the SINGLES dictionary.
    Pick k well-formed doubles whose both components are measured singles (comp_ok)."""
    idx = np.where(comp_ok)[0]
    rng = np.random.default_rng(0)
    pick = idx[np.argsort(-np.linalg.norm(Ed[idx], axis=1))][:k]   # strongest-effect doubles
    return [(int(i), doub[i]) for i in pick]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    os.makedirs(RESULTS, exist_ok=True)
    t0 = time.time()
    rows = []

    # ---------------- Norman (fast, rich) ----------------
    Es, Ed, sing, doub, comp_idx, comp_ok = load_norman()
    print(f"[norman] Es{Es.shape} Ed{Ed.shape}", flush=True)
    for di, dname in norman_targets(Es, Ed, sing, doub, comp_idx, comp_ok, k=6):
        d = Ed[di]
        for scheme in (None, SCHEME):
            r = calib_row(Es, d, dataset="Norman2019_K562_CRISPRa",
                          target=f"{dname} (held-out double)", weights_scheme=scheme,
                          extra=dict(target_kind="held-out double perturbation"))
            rows.append(r)
            print(f"  {dname:14s} {r['weighting']:18s} cos={r['reachable_cosine']:.3f} "
                  f"z={r['heldout_z']:.1f} DRF={r['dynamic_range_fraction']:.3f}", flush=True)
    pd.DataFrame(rows).to_csv(os.path.join(RESULTS, "deg_weighted_verdicts.csv"), index=False)
    print(f"[norman] done {time.time()-t0:.0f}s", flush=True)

    if args.quick:
        print("quick mode: skipping Tier-2", flush=True)
        _finalize(rows)
        return

    # ---------------- Tier-2 CD4 atlas (headline) ----------------
    import h5py
    def decode(a): return np.array([x.decode() if isinstance(x, bytes) else x for x in a])
    with h5py.File(os.path.join(DATA, "GWCD4i.DE_stats.h5ad"), "r") as h:
        var_gene = decode(h["var"]["gene_name"][:])
    gene_pos = {g: i for i, g in enumerate(var_gene)}
    pol = pd.read_csv(os.path.join(DATA, "Th2_Th1_polarization_signature_DE_results_full.suppl_table.csv"))
    piv = pol.pivot_table(index="variable", columns="contrast", values="zscore", aggfunc="mean").dropna()
    core = piv[np.sign(piv.iloc[:, 0]) == np.sign(piv.iloc[:, 1])]
    d_full = -1.0 * core.mean(axis=1)
    present = [g for g in d_full.index if g in gene_pos]
    d_t2 = d_full[present].values
    for cond in ("Rest", "Stim8hr", "Stim48hr"):
        E = np.load(os.path.join(CACHE, f"E_{cond}.npz"), allow_pickle=True)["E"].astype(np.float64)
        assert E.shape[1] == d_t2.shape[0], f"{cond} col mismatch"
        full = (cond == "Rest")     # full floor+ceiling+DRF calibration on the headline only
        print(f"[tier2] {cond} E{E.shape} (full_calibration={full})", flush=True)
        for scheme in (None, SCHEME):
            tt = time.time()
            if full:
                r = calib_row(E, d_t2, dataset="Tier2_CD4_CRISPRi",
                              target=f"toward_Th1 ({cond})", weights_scheme=scheme,
                              n_floor=N_FLOOR_T2, n_ceil=N_CEIL_T2, n_heldout=N_HELDOUT_T2,
                              extra=dict(target_kind="Th2->Th1 polarization signature"))
            else:
                r = robust_row(E, d_t2, dataset="Tier2_CD4_CRISPRi",
                               target=f"toward_Th1 ({cond})", weights_scheme=scheme,
                               n_heldout=N_HELDOUT_T2,
                               extra=dict(target_kind="Th2->Th1 polarization signature"))
            rows.append(r)
            drf = r["dynamic_range_fraction"]
            drf_s = f"{drf:.3f}" if drf == drf else "n/a"     # NaN-safe
            print(f"  {cond:9s} {r['weighting']:18s} cos={r['reachable_cosine']:.3f} "
                  f"z={r['heldout_z']:.1f} DRF={drf_s} ({time.time()-tt:.0f}s)", flush=True)
        pd.DataFrame(rows).to_csv(os.path.join(RESULTS, "deg_weighted_verdicts.csv"), index=False)
    _finalize(rows)
    print(f"ALL DONE {time.time()-t0:.0f}s", flush=True)


def _finalize(rows):
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(RESULTS, "deg_weighted_verdicts.csv"), index=False)
    # positive-control ceiling table (one row per dataset x weighting).
    # Use nan-aware first/mean so datasets where only the headline condition carries the full
    # calibration (Tier-2 Rest) still report their ceiling/floor.
    def _first_valid(s):
        v = s.dropna()
        return v.iloc[0] if len(v) else np.nan
    ceil = (df.groupby(["dataset", "weighting"])
              .agg(ceiling=("ceiling", _first_valid), ceiling_p25=("ceiling_p25", _first_valid),
                   ceiling_p75=("ceiling_p75", _first_valid), floor_mean=("floor_mean", "mean"))
              .reset_index())
    ceil.to_csv(os.path.join(RESULTS, "positive_control_ceiling.csv"), index=False)
    # calibration table (verdict placement)
    cal = df[["dataset", "target", "weighting", "reachable_cosine", "floor_mean",
              "ceiling", "dynamic_range_fraction", "heldout_z"]].copy()
    cal.to_csv(os.path.join(RESULTS, "calibration_dynamic_range.csv"), index=False)
    print(f"wrote {len(df)} rows -> deg_weighted_verdicts.csv, positive_control_ceiling.csv, "
          f"calibration_dynamic_range.csv", flush=True)


if __name__ == "__main__":
    main()
