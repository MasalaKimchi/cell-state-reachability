"""Atlas batch runner: 4 cell-state-transition targets x 3 conditions.
Full readout, full P. Signed decomposition (primary) + held-out-gene validation
(significance) + random-perturbation null (secondary). Writes each cell to
atlas_work/cell_<target>_<cond>.json as it completes so partial results survive.
Run from the workspace (where atlas_work/inputs.npz and reachability.py live).
"""
import numpy as np, json, os, sys, time
import os as _os
_REPO_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
_os.chdir(_REPO_ROOT); sys.path.insert(0, _REPO_ROOT)
import reachability as R

WORK = "analysis_cache/atlas_work"
z = np.load(f"{WORK}/inputs.npz", allow_pickle=True)
var_gene = z["var_gene"]
E = {c: z[f"E_{c}"].astype(np.float64) for c in ["Rest","Stim8hr","Stim48hr"]}
GENE = {c: z[f"gene_{c}"] for c in ["Rest","Stim8hr","Stim48hr"]}
TARGETS = {n: z[f"t_{n}"] for n in ["toward_Th1","toward_Th2","toward_younger","toward_older"]}

# Null counts. Defaults (15, 8) are UNCHANGED so a rerun writes byte-identical JSON to
# before this change; only the SOLVE MECHANISM is faster now (Gram reuse + process
# parallelism in reachability.py), not the statistics. Because solves are ~6x cheaper you
# can now AFFORD a stronger held-out null via env override (this DOES change the numbers):
#   N_HELDOUT=60 python scripts/run_atlas.py      # stronger null (different z, by design)
N_HELDOUT = int(_os.environ.get("N_HELDOUT", "15"))   # shuffled-target held-out null count
N_PERT    = int(_os.environ.get("N_PERT", "8"))       # random-perturbation null count
N_JOBS    = int(_os.environ.get("REACH_N_JOBS", "-1"))  # -1 = all cores; 1 = serial (output-identical)

def point_estimate(name, d_full, cond):
    """Fast hero verdict: signed decomposition only (2 NNLS)."""
    Em = E[cond]; idx = np.where(d_full != 0)[0]
    mask = np.zeros(len(var_gene), bool); mask[idx] = True
    d = d_full.copy(); t0 = time.time()
    s = R.signed_reachability(Em, d, hvg_mask=mask)
    # canonical residual_norm = ||d_masked - fitted_lof|| / ||d_masked||  (matches reachability())
    dm = d[mask]; dmn = np.linalg.norm(dm)
    residual_norm = float(np.linalg.norm(dm - s.fitted_lof) / dmn) if dmn else 0.0
    gene_axis = GENE[cond]   # perturbed-gene name per generator row
    # top greedy nominations (fast OMP path) for the recipe
    try:
        spec = R.reachability_spectrum(Em, d, k_max=10, hvg_mask=mask, refit_full=False)
        greedy = [str(gene_axis[int(g)]) for g in spec["order"][:10]]
    except Exception as e:
        greedy = [f"<greedy failed: {e}>"]
    return dict(target=name, condition=cond, n_readout=int(mask.sum()), n_generators=int(Em.shape[0]),
        reachable_cosine=s.lof_cosine, residual_norm=residual_norm,
        lof_fraction=s.lof_fraction, gof_fraction=s.gof_fraction,
        neither_fraction=s.neither_fraction, signed_cosine=s.signed_cosine, kkt_viol=s.cert_max_violation,
        lof_support_size=int(s.lof_support.size), gof_support_size=int(s.gof_support.size),
        greedy_order=greedy, seconds=round(time.time()-t0,1))

def nulls(name, d_full, cond, seed=0, lof_cosine=None):
    """Significance: held-out-gene validation + random-perturbation null.

    lof_cosine : the LOF reachable cosine already computed by point_estimate (persisted as
        `reachable_cosine`). Passed in to avoid recomputing an identical ~46s signed_reachability
        just to read s.lof_cosine back for the perturbation-null z. Falls back to recomputing
        only if not supplied (e.g. nulls() called standalone).
    """
    Em = E[cond]; idx = np.where(d_full != 0)[0]
    mask = np.zeros(len(var_gene), bool); mask[idx] = True
    d = d_full.copy(); t0 = time.time()
    ho = R.held_out_gene_validation(Em, d, hvg_mask=mask, n_shuffles=N_HELDOUT, seed=seed,
                                    n_jobs=N_JOBS)
    rng = np.random.default_rng(seed + 7)
    Asub = Em[:, idx].T; dsub = d[idx]; dn = np.linalg.norm(dsub); K, P = Asub.shape
    from scipy.optimize import nnls
    if lof_cosine is None:                        # standalone fallback: recompute
        lof_cosine = R.signed_reachability(Em, d, hvg_mask=mask).lof_cosine
    pert = []
    for _ in range(N_PERT):
        perm = np.argsort(rng.random((K, P)), axis=0)
        Ap = np.take_along_axis(Asub, perm, axis=0)
        w, _ = nnls(Ap, dsub); fit = Ap @ w
        pert.append(float(fit @ dsub / (np.linalg.norm(fit) * dn + 1e-12)))
    pert = np.array(pert)
    return dict(heldout_cosine=ho.held_out_cosine, heldout_null_mean=ho.null_mean,
        heldout_null_std=ho.null_std, heldout_z=ho.z,
        pert_null_mean=float(pert.mean()), pert_null_std=float(pert.std()),
        pert_null_z=float((lof_cosine - pert.mean())/(pert.std()+1e-9)),
        null_seconds=round(time.time()-t0,1))

def main():
    conds = ["Rest","Stim8hr","Stim48hr"]
    order = [(t,c) for t in TARGETS for c in conds]
    # PASS 1: point estimates (hero verdicts) — all 12 land fast
    print("=== PASS 1: point estimates ===", flush=True)
    for i,(name,cond) in enumerate(order,1):
        fp=f"{WORK}/point_{name}_{cond}.json"
        if os.path.exists(fp): print(f"skip point {name} {cond}",flush=True); continue
        out=point_estimate(name,TARGETS[name],cond); json.dump(out,open(fp,"w"),indent=2)
        print(f"[P{i}/12] {name} {cond}: reach={out['reachable_cosine']:.3f} LOF={out['lof_fraction']:.2f} "
              f"GOF={out['gof_fraction']:.2f} neither={out['neither_fraction']:.2f} greedy={out['greedy_order'][:5]} ({out['seconds']}s)",flush=True)
    print("=== PASS 1 COMPLETE ===", flush=True)
    # PASS 2: nulls — augment each cell
    print("=== PASS 2: nulls ===", flush=True)
    for i,(name,cond) in enumerate(order,1):
        fp=f"{WORK}/cell_{name}_{cond}.json"
        if os.path.exists(fp): print(f"skip null {name} {cond}",flush=True); continue
        pt=json.load(open(f"{WORK}/point_{name}_{cond}.json"))
        nl=nulls(name,TARGETS[name],cond,lof_cosine=pt["reachable_cosine"])  # reuse PASS-1 cosine
        pt.update(nl); json.dump(pt,open(fp,"w"),indent=2)
        print(f"[N{i}/12] {name} {cond}: heldout_z={pt['heldout_z']:.1f} pert_z={pt['pert_null_z']:.1f} ({pt['null_seconds']}s)",flush=True)
    print("ATLAS COMPLETE", flush=True)

if __name__ == "__main__":
    main()
