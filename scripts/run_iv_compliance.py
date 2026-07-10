"""Instrumental-variables / compliance robustness of the reachability verdict.

Reframes the CRISPRi Perturb-seq oracle in causal-inference terms:
  * instrument Z  = guide assignment (randomized in the pool)
  * treatment  T  = realized on-target knockdown of the perturbed gene (imperfect,
                    heterogeneous compliance; measured in guide_kd_efficiency)
  * outcome    Y  = transcriptomic shift  => the effect vectors E[g,:] are INTENT-TO-TREAT
                    (the effect of *assigning* guide g, not of achieving knockdown)

Three verdicts per (target, condition), all via signed_reachability on the SAME readout:
  A. ITT           — the full dictionary (as published; current headline).
  B. valid-instr.  — restrict to perturbations with >=1 significant on-target knockdown
                     (exclusion-restriction robustness: drop generators that move the
                     transcriptome without demonstrably hitting their intended target).
  Recipe/LATE and compliance-headroom are derived downstream from B's weights + compliance.

Writes atlas_work/iv_cell_<target>_<cond>.json as each lands (partials survive).
Run from the workspace (needs atlas_work/inputs.npz, atlas_work/first_stage_compliance.csv,
reachability.py).
"""
import numpy as np, pandas as pd, json, os, sys, time
import os as _os
_REPO_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
_os.chdir(_REPO_ROOT); sys.path.insert(0, _REPO_ROOT)
import reachability as R

WORK = "analysis_cache/atlas_work"
z = np.load(f"{WORK}/inputs.npz", allow_pickle=True)
var_gene = z["var_gene"]
CONDS = ["Rest", "Stim8hr", "Stim48hr"]
E     = {c: z[f"E_{c}"].astype(np.float64) for c in CONDS}
GENE  = {c: z[f"gene_{c}"] for c in CONDS}
TARGETS = {n: z[f"t_{n}"] for n in ["toward_Th1", "toward_Th2", "toward_younger", "toward_older"]}

fs = pd.read_csv(f"{WORK}/first_stage_compliance.csv")
# per (cond) lookup: gene -> (valid_instrument, pi_cellwt)
VALID = {}; PI = {}
for c in CONDS:
    sub = fs[fs["cond"] == c].set_index("gene")
    VALID[c] = sub["valid_instrument"].to_dict()
    PI[c]    = sub["pi_cellwt"].to_dict()


def solve(name, d_full, cond, restrict):
    """One signed_reachability solve; restrict=True keeps only valid instruments."""
    Em = E[cond]; gene_axis = GENE[cond]
    idx = np.where(d_full != 0)[0]
    mask = np.zeros(len(var_gene), bool); mask[idx] = True
    d = d_full.copy()
    if restrict:
        keep = np.array([bool(VALID[cond].get(str(g), False)) for g in gene_axis])
    else:
        keep = np.ones(len(gene_axis), bool)
    Ek = Em[keep]; gk = gene_axis[keep]
    s = R.signed_reachability(Ek, d, hvg_mask=mask)
    dm = d[mask]; dmn = np.linalg.norm(dm)
    residual_norm = float(np.linalg.norm(dm - s.fitted_lof) / dmn) if dmn else 0.0
    # greedy recipe (order + refit weights) on the restricted dictionary
    spec = R.reachability_spectrum(Ek, d, k_max=12, hvg_mask=mask, refit_full=True)
    order = [int(g) for g in spec["order"][:12]]
    recipe = [str(gk[g]) for g in order]
    recipe_pi = [float(PI[cond].get(str(gk[g]), np.nan)) for g in order]
    return dict(
        target=name, condition=cond, arm=("valid_instrument" if restrict else "ITT"),
        n_readout=int(mask.sum()), n_generators=int(Ek.shape[0]),
        n_generators_full=int(Em.shape[0]), n_dropped=int((~keep).sum()),
        reachable_cosine=float(s.lof_cosine), residual_norm=residual_norm,
        lof_fraction=float(s.lof_fraction), gof_fraction=float(s.gof_fraction),
        neither_fraction=float(s.neither_fraction), signed_cosine=float(s.signed_cosine),
        kkt_viol=float(s.cert_max_violation),
        recipe=recipe, recipe_compliance=recipe_pi,
    )


def late_rescaling_invariance(name, d_full, cond):
    """Convex-cone invariance check: rescale each generator row by 1/pi (the ITT->LATE
    Wald/compliance rescaling) and confirm the reachable cosine is unchanged. A convex
    cone is invariant to positive per-generator rescaling, so this must return ~0 to
    machine precision (it relabels the recipe weights, not the verdict)."""
    Em = E[cond]; gene_axis = GENE[cond]
    idx = np.where(d_full != 0)[0]
    mask = np.zeros(len(var_gene), bool); mask[idx] = True
    d = d_full.copy()
    pi = np.array([PI[cond].get(str(g), np.nan) for g in gene_axis])
    ok = np.isfinite(pi) & (pi > 0.05)          # valid, well-measured instruments
    scale = np.where(ok, 1.0 / np.where(ok, pi, 1.0), 1.0)   # 1/pi, floored for weak/invalid
    cos_itt = float(R.reachability(Em, d, hvg_mask=mask).reachable_cosine)
    cos_late = float(R.reachability((Em.T * scale).T, d, hvg_mask=mask).reachable_cosine)
    return dict(cos_ITT=cos_itt, cos_LATE=cos_late, abs_delta=abs(cos_itt - cos_late),
                scale_median=float(np.median(scale)), scale_max=float(scale.max()))


def write_invariance_csv(path=f"{WORK}/late_rescaling_invariance.csv"):
    """Compute the ITT->LATE rescaling invariance for all 12 (target, condition) cells
    via late_rescaling_invariance() and write the per-cell table to CSV. This is the
    routine that produces atlas_work/late_rescaling_invariance.csv cited in
    CAUSAL.md (IV/compliance layer) - self-contained (no dependence on the per-cell JSON pipeline)."""
    rows = []
    for name in TARGETS:
        for cond in CONDS:
            lr = late_rescaling_invariance(name, TARGETS[name], cond)
            rows.append(dict(target=name.replace("toward_", ""), cond=cond,
                             cos_ITT=lr["cos_ITT"], cos_LATE=lr["cos_LATE"],
                             abs_delta=lr["abs_delta"], scale_median=lr["scale_median"],
                             scale_max=lr["scale_max"]))
    df = pd.DataFrame(rows)
    df.to_csv(path, index=False)
    print(f"wrote {path}: max |Δcosine| = {df.abs_delta.max():.2e} over {len(df)} cells",
          flush=True)
    return df


def main():
    order = [(t, c) for t in TARGETS for c in CONDS]
    for i, (name, cond) in enumerate(order, 1):
        fp = f"{WORK}/iv_cell_{name}_{cond}.json"
        if os.path.exists(fp):
            print(f"skip {name} {cond}", flush=True); continue
        t0 = time.time()
        itt = solve(name, TARGETS[name], cond, restrict=False)
        val = solve(name, TARGETS[name], cond, restrict=True)
        late = late_rescaling_invariance(name, TARGETS[name], cond)   # ITT->LATE invariance
        out = {"ITT": itt, "valid_instrument": val, "late_rescaling": late,
               "delta_cosine": val["reachable_cosine"] - itt["reachable_cosine"],
               "seconds": round(time.time() - t0, 1)}
        json.dump(out, open(fp, "w"), indent=2)
        print(f"[{i}/12] {name} {cond}: ITT cos={itt['reachable_cosine']:.3f} "
              f"-> valid cos={val['reachable_cosine']:.3f} "
              f"(drop-Δ={out['delta_cosine']:+.4f}, dropped {val['n_dropped']}/{itt['n_generators_full']}; "
              f"LATE-rescale |Δ|={late['abs_delta']:.1e}) "
              f"({out['seconds']}s)", flush=True)
    print("IV COMPLIANCE COMPLETE", flush=True)


if __name__ == "__main__":
    if "--invariance-only" in sys.argv:
        write_invariance_csv()          # just the LATE-rescaling CSV (fast: 12 x 2 solves)
    else:
        main()                          # full per-cell JSON pipeline (ITT + valid + LATE)
        write_invariance_csv()          # and the invariance summary CSV
