"""Pass 2 (lean): held-out-gene validation z per atlas cell — the primary
significance test. Reads atlas_work/point_*.json, augments each into
atlas_work/cell_<target>_<cond>.json. Restartable (skips existing cell_*).
The expensive random-perturbation null is validated separately on the headline
(pert_z=-18.4) and omitted here for tractability.
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
TARGETS = {n: z[f"t_{n}"] for n in ["toward_Th1","toward_Th2","toward_younger","toward_older"]}
N_HELDOUT = 8   # each shuffle is a full-readout NNLS refit (~14s); z stays >>3 (null is tight)

def held(name, d_full, cond, seed=0):
    Em = E[cond]; idx = np.where(d_full != 0)[0]
    mask = np.zeros(len(var_gene), bool); mask[idx] = True
    t0 = time.time()
    ho = R.held_out_gene_validation(Em, d_full.copy(), hvg_mask=mask, n_shuffles=N_HELDOUT, seed=seed)
    return dict(heldout_cosine=float(ho.held_out_cosine), heldout_null_mean=float(ho.null_mean),
                heldout_null_std=float(ho.null_std), heldout_z=float(ho.z),
                n_heldout_shuffles=N_HELDOUT, null_seconds=round(time.time()-t0,1))

def main():
    conds=["Rest","Stim8hr","Stim48hr"]
    order=[(t,c) for t in TARGETS for c in conds]
    for i,(name,cond) in enumerate(order,1):
        fp=f"{WORK}/cell_{name}_{cond}.json"
        if os.path.exists(fp): print(f"skip {name} {cond}",flush=True); continue
        pt=json.load(open(f"{WORK}/point_{name}_{cond}.json"))
        pt.update(held(name,TARGETS[name],cond))
        json.dump(pt,open(fp,"w"),indent=2)
        print(f"[N{i}/12] {name} {cond}: heldout_cos={pt['heldout_cosine']:.3f} "
              f"heldout_z={pt['heldout_z']:.1f} ({pt['null_seconds']}s)",flush=True)
    print("NULLS COMPLETE",flush=True)

if __name__=="__main__":
    main()
