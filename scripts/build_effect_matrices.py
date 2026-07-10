"""Build row/column-aligned effect matrices E_K562 and E_RPE1 from the CZI
harmonized Replogle-2022 essential-gene Perturb-seq files.

Effect definition (per perturbation g, per gene):
    E[g, :] = mean(X | condition==g) - mean(X | control cells)
where X is the LOG-NORMALIZED expression already stored in the files (verified:
non-negative, log-compressed 0.24-6.3, ~2091 row sums; NOT raw counts, so we do
NOT re-normalize). This is a pseudobulk log fold-change vs the non-targeting
control pool, computed IDENTICALLY in both cell types. No cross-perturbation
per-gene centering (that would inject the rank-deficiency identity we flagged in
the Norman analysis); each perturbation's effect vector is independent.

Group means are computed with a single sparse matmul: an (n_groups x n_cells)
row-normalized indicator G times X gives all per-condition means at once.

Alignment:
  * genes  -> intersection of var.gene_id (Ensembl, stable across cell types)
  * perturbations -> single-gene targets present in BOTH cell types
Both matrices come out (P_shared, G_shared) with identical row and column order,
so E_K562 and E_RPE1 are directly comparable element-wise.
"""
import numpy as np, scipy.sparse as sp, anndata as ad, json, sys, os
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def load_effect(path, tag):
    print(f"[{tag}] loading {path}", flush=True)
    A = ad.read_h5ad(path)
    print(f"[{tag}] {A.shape[0]:,} cells x {A.shape[1]:,} genes", flush=True)
    X = A.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    X = X.astype(np.float64).tocsr()          # already log-normalized; use as-is
    cond = A.obs["condition"].astype(str).values
    ctrl = np.asarray(A.obs["control"]).astype(int).ravel()
    gene_id = np.asarray(A.var["gene_id"]).astype(str) if "gene_id" in A.var else np.asarray(A.var_names).astype(str)
    gene_name = np.asarray(A.var["feature_name"]).astype(str) if "feature_name" in A.var else gene_id.copy()

    # control pool mean
    ctrl_mask = ctrl == 1
    print(f"[{tag}] control cells: {ctrl_mask.sum():,}", flush=True)
    mu_ctrl = np.asarray(X[ctrl_mask].mean(axis=0)).ravel()

    # single-gene targets only: 'GENE+ctrl'
    singles = [c for c in np.unique(cond) if c.endswith("+ctrl") and c[:-5] not in ("", "ctrl")]
    genes = [c[:-5] for c in singles]
    # build (n_groups x n_cells) row-normalized indicator, one matmul for all means
    col = {c: j for j, c in enumerate(singles)}
    rows, cols = [], []
    for i, c in enumerate(cond):
        j = col.get(c)
        if j is not None:
            rows.append(j); cols.append(i)
    counts = np.bincount(rows, minlength=len(singles)).astype(float)
    data = np.ones(len(rows))
    Gind = sp.csr_matrix((data, (rows, cols)), shape=(len(singles), X.shape[0]))
    means = np.asarray((Gind @ X).todense()) / counts[:, None]   # (n_groups, n_genes)
    eff = {genes[i]: means[i] - mu_ctrl for i in range(len(singles))}
    ncells = {genes[i]: int(counts[i]) for i in range(len(singles))}
    print(f"[{tag}] built effect vectors for {len(eff):,} single-gene perturbations", flush=True)
    return dict(eff=eff, ncells=ncells, gene_id=gene_id, gene_name=gene_name)

def main():
    K = load_effect("analysis_cache/czi_data/k562_essential.h5ad", "K562")
    R = load_effect("analysis_cache/czi_data/rpe1_essential.h5ad", "RPE1")

    # shared genes by Ensembl id
    gk = {gid: i for i, gid in enumerate(K["gene_id"])}
    gr = {gid: i for i, gid in enumerate(R["gene_id"])}
    shared_genes = [gid for gid in K["gene_id"] if gid in gr]
    ik = [gk[g] for g in shared_genes]
    ir = [gr[g] for g in shared_genes]
    gene_name_shared = [K["gene_name"][gk[g]] for g in shared_genes]
    print(f"shared genes: {len(shared_genes):,}", flush=True)

    # shared perturbations
    shared_perts = sorted(set(K["eff"]) & set(R["eff"]))
    print(f"shared perturbations: {len(shared_perts):,}", flush=True)

    Ek = np.vstack([K["eff"][g][ik] for g in shared_perts])
    Er = np.vstack([R["eff"][g][ir] for g in shared_perts])
    nk = np.array([K["ncells"][g] for g in shared_perts])
    nr = np.array([R["ncells"][g] for g in shared_perts])
    print(f"E_K562 {Ek.shape}  E_RPE1 {Er.shape}", flush=True)

    np.savez_compressed(
        "analysis_cache/czi_data/cross_celltype_effects.npz",
        E_K562=Ek.astype(np.float32), E_RPE1=Er.astype(np.float32),
        perturbations=np.array(shared_perts),
        gene_id=np.array(shared_genes), gene_name=np.array(gene_name_shared),
        n_cells_K562=nk, n_cells_RPE1=nr,
    )
    meta = dict(n_perturbations=len(shared_perts), n_genes=len(shared_genes),
                E_K562_shape=list(Ek.shape), E_RPE1_shape=list(Er.shape),
                effect="pseudobulk log1p(CP10k) mean(pert) - mean(control)",
                source="s3://cz-benchmarks-data Replogle 2022 K562/RPE1 essential")
    json.dump(meta, open("analysis_cache/czi_data/cross_celltype_effects_meta.json", "w"), indent=2)
    print("SAVED analysis_cache/czi_data/cross_celltype_effects.npz", flush=True)
    print(json.dumps(meta, indent=2))

if __name__ == "__main__":
    main()
