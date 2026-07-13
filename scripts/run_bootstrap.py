"""Robustness layer: gene-panel SUBSAMPLING stability on the reachable cosine +
signed decomposition. Draw random 85%-without-replacement subsets of the readout
genes, refit the signed decomposition, and report the 2.5-97.5 percentile spread.
This is the HONEST substitute for leave-one-donor-out: the effect matrix is
donor-collapsed (only cross-donor QC correlations survive; no per-donor effect
vectors are local), so true LODO is impossible. The gene-panel subsample instead
asks how sensitive the LOF/GOF/neither verdict is to the specific signature gene
set. (Subsampling without replacement, not with-replacement bootstrap, so NNLS
stays fast — duplicated columns triple NNLS time.)
Runs the 3 Rest 'PARTIALLY REACHABLE' headline cells.
"""
import numpy as np, json, sys, time
import os as _os
_REPO_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
_os.chdir(_REPO_ROOT); sys.path.insert(0, _REPO_ROOT)
import reachability as R
from scipy.optimize import nnls

z=np.load("analysis_cache/atlas_work/inputs.npz",allow_pickle=True)
var_gene=z["var_gene"]
N_SUB=int(_os.environ.get("N_SUB","20")); FRAC=0.85
# Process-parallelism over the independent subsamples (-1 = all cores, 1 = serial). Each
# iteration subsamples the GENE axis, so the fit matrix A=Asub[s] changes every iteration --
# Gram-reuse does NOT apply here (the Gram matrix is not fixed). Only the n_jobs half of the
# fast path is valid, and it is output-preserving: all subsample index sets are drawn up
# front, so results are independent of worker count. Override with REACH_N_JOBS=1 for serial.
N_JOBS=int(_os.environ.get("REACH_N_JOBS","-1"))

# Fork-worker state: the (fixed) full sub-dictionary lives in a module global so children
# inherit it copy-on-write rather than re-pickling it per task.
_BW={}

def _bw_init(Asub, dsub):
    try:
        from threadpoolctl import threadpool_limits; threadpool_limits(1)
    except Exception:
        pass
    _BW["Asub"]=Asub; _BW["dsub"]=dsub

def _bw_solve(s):
    """One 85%-gene-subsample signed refit; s = the selected gene indices."""
    A=_BW["Asub"][s]; dd=_BW["dsub"][s]; dn=np.linalg.norm(dd); tot=dn*dn
    wl,_=nnls(A,dd); fl=A@wl; r=dd-fl
    wg,_=nnls(-A,r); fg=-A@wg; rho=r-fg
    return (float(fl@dd/(np.linalg.norm(fl)*dn+1e-12)),
            float((tot-np.dot(r,r))/tot), float((np.dot(r,r)-np.dot(rho,rho))/tot),
            float(np.dot(rho,rho)/tot))

def _bw_map(subsets):
    """Map _bw_solve over subsample index sets; parallel via fork pool, serial fallback."""
    n=R._resolve_n_jobs(N_JOBS)
    if n<=1 or len(subsets)<=1:
        return [_bw_solve(s) for s in subsets]
    import multiprocessing as mp
    try:
        ctx=mp.get_context("fork")
        with ctx.Pool(n, initializer=_bw_init, initargs=(_BW["Asub"],_BW["dsub"])) as pool:
            return pool.map(_bw_solve, subsets)
    except (OSError, ValueError, ImportError):
        return [_bw_solve(s) for s in subsets]

def sub_cell(tkey, ckey, seed=0):
    E=z[f"E_{ckey}"].astype(np.float64); d=z[f"t_{tkey}"]
    idx=np.where(d!=0)[0]; mask=np.zeros(len(var_gene),bool); mask[idx]=True
    s0=R.signed_reachability(E,d.copy(),hvg_mask=mask)          # point estimate
    Asub=E[:,idx].T; dsub=d[idx]; K=len(idx); rng=np.random.default_rng(seed)
    m=int(round(FRAC*K))
    # Draw all subsample sets up front so output is identical regardless of worker count.
    subsets=[rng.choice(K,m,replace=False) for _ in range(N_SUB)]
    _bw_init(Asub, dsub)                                        # also primes serial path
    t0=time.time()
    results=_bw_map(subsets)
    rc=np.array([x[0] for x in results]); lof=np.array([x[1] for x in results])
    gof=np.array([x[2] for x in results]); nei=np.array([x[3] for x in results])
    def ci(a): return [round(float(np.percentile(a,2.5)),4),round(float(np.percentile(a,97.5)),4)]
    return dict(target=tkey,condition=ckey,method=f"{int(FRAC*100)}% gene subsample x{N_SUB}",
        point_reach_cos=round(float(s0.lof_cosine),4),point_lof=round(float(s0.lof_fraction),4),
        point_gof=round(float(s0.gof_fraction),4),point_neither=round(float(s0.neither_fraction),4),
        sub_reach_cos_mean=round(float(rc.mean()),4),sub_reach_cos_ci=ci(rc),
        sub_lof_mean=round(float(lof.mean()),4),sub_lof_ci=ci(lof),
        sub_gof_ci=ci(gof),sub_neither_ci=ci(nei),
        n_sub=N_SUB,seconds=round(time.time()-t0,0))

if __name__=="__main__":
    res=[]
    for tk in ["toward_Th1","toward_Th2","toward_younger"]:
        r=sub_cell(tk,"Rest",seed=0); res.append(r)
        print(f"{tk} Rest: reach_cos={r['point_reach_cos']} CI={r['sub_reach_cos_ci']} "
              f"| LOF={r['point_lof']} CI={r['sub_lof_ci']} ({r['seconds']}s)",flush=True)
        json.dump(res,open("analysis_cache/atlas_work/bootstrap_ci.json","w"),indent=2)
    print("BOOTSTRAP COMPLETE",flush=True)
