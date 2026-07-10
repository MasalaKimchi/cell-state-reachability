import sys, json, numpy as np
sys.path.insert(0, "..")                     # project root: reachability.py
import reachability as rr

# Effect-matrix checkpoint (built by build_effect_matrices.py from the two CZI h5ad files)
d = np.load("../czi_data/cross_celltype_effects.npz", allow_pickle=True)
E_K562, E_RPE1 = d["E_K562"].astype(float), d["E_RPE1"].astype(float)
perts   = d["perturbations"].astype(str)
gene_id = d["gene_id"].astype(str)
gene_nm = d["gene_name"].astype(str)
P, G = E_K562.shape
print(f"E_K562 {E_K562.shape}  E_RPE1 {E_RPE1.shape}")
print(f"{P} shared single-gene perturbations x {G} shared readout genes")

name2col = {g:i for i,g in enumerate(gene_nm)}
on_k, on_r = [], []
for pi, p in enumerate(perts):
    j = name2col.get(p)
    if j is not None:
        on_k.append(E_K562[pi, j]); on_r.append(E_RPE1[pi, j])
on_k, on_r = np.array(on_k), np.array(on_r)
print(f"on-target self-effect checked on {len(on_k)} perturbations")
print(f"  K562: mean {on_k.mean():+.3f}, fraction negative {np.mean(on_k<0):.1%}")
print(f"  RPE1: mean {on_r.mean():+.3f}, fraction negative {np.mean(on_r<0):.1%}")

def cos(a, b):
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    return float(a @ b / (na * nb)) if na and nb else 0.0

rng = np.random.default_rng(0)
matched   = np.array([cos(E_K562[i], E_RPE1[i]) for i in range(P)])
null_shuf = np.array([cos(E_K562[i], E_RPE1[i][rng.permutation(G)]) for i in range(P)])
jj = rng.permutation(P); jj[jj==np.arange(P)] = (jj[jj==np.arange(P)]+1) % P
null_mm   = np.array([cos(E_K562[i], E_RPE1[jj[i]]) for i in range(P)])

print(f"matched cosine        median {np.median(matched):+.3f}")
print(f"shuffled-gene null    median {np.median(null_shuf):+.3f}  (95th {np.percentile(null_shuf,95):+.3f})")
print(f"mismatched-pert null  median {np.median(null_mm):+.3f}  (95th {np.percentile(null_mm,95):+.3f})")

# (a) per-perturbation gene-specificity z
K = 200
spec_z = np.zeros(P); spec_p = np.zeros(P)
rng = np.random.default_rng(1)
for i in range(P):
    others = rng.choice(np.delete(np.arange(P), i), size=K, replace=False)
    mm = np.array([cos(E_K562[i], E_RPE1[o]) for o in others])
    spec_z[i] = (matched[i]-mm.mean())/mm.std() if mm.std() else 0.0
    spec_p[i] = (mm >= matched[i]).mean()

# (b) deflate common essential-stress direction
def deflate(E, c):
    cu = c/np.linalg.norm(c); return E - np.outer(E@cu, cu)
common_cos = cos(E_K562.mean(0), E_RPE1.mean(0))
Ek_d, Er_d = deflate(E_K562, E_K562.mean(0)), deflate(E_RPE1, E_RPE1.mean(0))
matched_defl = np.array([cos(Ek_d[i], Er_d[i]) for i in range(P)])

print(f"common essential-stress direction cosine: {common_cos:+.3f}")
print(f"gene-specific at p<0.05: {np.mean(spec_p<0.05):.1%};  z>2: {np.mean(spec_z>2):.1%}")
print(f"matched cosine after deflation: median {np.median(matched_defl):+.3f}; still positive for {np.mean(matched_defl>0):.1%}")

idx = np.arange(P)
def fit(E_basis, dvec):
    r = rr.reachability(E_basis, dvec)
    return r.residual_norm, r.reachable_cosine

win_k_res=np.zeros(P); win_r_res=np.zeros(P); crs_kr_res=np.zeros(P); crs_rk_res=np.zeros(P)
win_k_cos=np.zeros(P); win_r_cos=np.zeros(P); crs_kr_cos=np.zeros(P); crs_rk_cos=np.zeros(P)
for i in idx:
    m = idx != i
    win_k_res[i], win_k_cos[i] = fit(E_K562[m], E_K562[i])   # K562 target, K562 basis
    win_r_res[i], win_r_cos[i] = fit(E_RPE1[m], E_RPE1[i])   # RPE1 target, RPE1 basis
    crs_kr_res[i], crs_kr_cos[i] = fit(E_RPE1[m], E_K562[i]) # K562 target, RPE1 basis
    crs_rk_res[i], crs_rk_cos[i] = fit(E_K562[m], E_RPE1[i]) # RPE1 target, K562 basis

print("median cone residual (lower = more reachable):")
print(f"  within K562 {np.median(win_k_res):.3f}   cross K562->RPE1 basis {np.median(crs_kr_res):.3f}")
print(f"  within RPE1 {np.median(win_r_res):.3f}   cross RPE1->K562 basis {np.median(crs_rk_res):.3f}")

def null_band(E_basis, n=60, seed=0):
    rng = np.random.default_rng(seed); Pn = E_basis.shape[0]
    sel = rng.choice(Pn, n, replace=False); out=[]
    for i in sel:
        dp = E_basis[i][rng.permutation(E_basis.shape[1])]
        out.append(rr.reachability(np.delete(E_basis,i,0), dp).reachable_cosine)
    return np.array(out)

nb_k, nb_r = null_band(E_K562,60,0), null_band(E_RPE1,60,1)
thr_k, thr_r = np.percentile(nb_k,95), np.percentile(nb_r,95)
v_win_k, v_crs_kr = win_k_cos>thr_k, crs_kr_cos>thr_r   # K562 targets: within(K basis), cross(R basis)
v_win_r, v_crs_rk = win_r_cos>thr_r, crs_rk_cos>thr_k   # RPE1 targets
print(f"null 95th cosine: K562 basis {thr_k:.3f}, RPE1 basis {thr_r:.3f}")
print(f"reachable within/cross:  K562 tgt {v_win_k.mean():.0%}/{v_crs_kr.mean():.0%}   RPE1 tgt {v_win_r.mean():.0%}/{v_crs_rk.mean():.0%}")
print(f"verdict agreement:  K562 tgt {(v_win_k==v_crs_kr).mean():.1%}   RPE1 tgt {(v_win_r==v_crs_rk).mean():.1%}")

def support_global(E_basis, i, target_vec, k=10):
    m = idx != i; gg = idx[m]
    return set(gg[rr.reachability(E_basis[m], target_vec).support[:k]].tolist())

def jac(a,b): return len(a&b)/len(a|b) if (a|b) else 0.0
jac_q1=np.zeros(P); jac_q2k=np.zeros(P); jac_q2r=np.zeros(P)
for i in idx:
    sk = support_global(E_K562, i, E_K562[i])   # K562 target, native
    sr = support_global(E_RPE1, i, E_RPE1[i])   # RPE1 target, native
    skx = support_global(E_RPE1, i, E_K562[i])  # K562 target from RPE1 basis
    srx = support_global(E_K562, i, E_RPE1[i])  # RPE1 target from K562 basis
    jac_q1[i]=jac(sk,sr); jac_q2k[i]=jac(sk,skx); jac_q2r[i]=jac(sr,srx)

rng=np.random.default_rng(0); nn=[]
for _ in range(20000):
    a=set(rng.choice(P-1,10,replace=False)); b=set(rng.choice(P-1,10,replace=False)); nn.append(jac(a,b))
null95=np.percentile(nn,95)
print(f"random-null Jaccard: mean {np.mean(nn):.3f}, 95th {null95:.3f}")
print(f"same-target recipe overlap: median {np.median(jac_q1):.3f}  ({np.mean(jac_q1>null95):.0%} above null-95th)")
print(f"cross-basis recipe overlap: median {np.median(jac_q2k):.3f} (K562 tgt) / {np.median(jac_q2r):.3f} (RPE1 tgt)")