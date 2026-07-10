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
N_SUB=20; FRAC=0.85

def sub_cell(tkey, ckey, seed=0):
    E=z[f"E_{ckey}"].astype(np.float64); d=z[f"t_{tkey}"]
    idx=np.where(d!=0)[0]; mask=np.zeros(len(var_gene),bool); mask[idx]=True
    s0=R.signed_reachability(E,d.copy(),hvg_mask=mask)          # point estimate
    Asub=E[:,idx].T; dsub=d[idx]; K=len(idx); rng=np.random.default_rng(seed)
    m=int(round(FRAC*K))
    rc=[]; lof=[]; gof=[]; nei=[]
    t0=time.time()
    for b in range(N_SUB):
        s=rng.choice(K,m,replace=False); A=Asub[s]; dd=dsub[s]; dn=np.linalg.norm(dd); tot=dn*dn
        wl,_=nnls(A,dd); fl=A@wl; r=dd-fl
        wg,_=nnls(-A,r); fg=-A@wg; rho=r-fg
        rc.append(float(fl@dd/(np.linalg.norm(fl)*dn+1e-12)))
        lof.append(float((tot-np.dot(r,r))/tot)); gof.append(float((np.dot(r,r)-np.dot(rho,rho))/tot))
        nei.append(float(np.dot(rho,rho)/tot))
    rc=np.array(rc); lof=np.array(lof); gof=np.array(gof); nei=np.array(nei)
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
