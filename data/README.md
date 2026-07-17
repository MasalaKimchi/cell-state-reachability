# External data

Scientific data are intentionally local and ignored by Git. The maintained repository
ships frozen selected evidence, not the 16.8 GB source matrix.

## Primary source

Zhu et al. 2025, genome-scale CRISPRi Perturb-seq in primary human CD4 cells:

- doi:10.64898/2025.12.23.696273
- SRA SRP643211 / GEO GSE314342

- `GWCD4i.DE_stats.h5ad`
- 33,983 perturbation-condition profiles × 10,282 readout genes
- relevant layers include `log_fc`, `zscore`, `p_value`, `adj_p_value`, `baseMean`, and
  `lfcSE`

Obtain the artifact through the CZI Virtual Cells Platform:

```bash
vcp data search "Primary Human CD4+ T Cell Perturb-seq" --exact
```

Registration and CLI documentation:
https://chanzuckerberg.github.io/vcp-cli/usage/data.html

Open supplementary tables are also available from the source-analysis repository:
https://github.com/emdann/GWT_perturbseq_analysis_2025/tree/master/metadata/suppl_tables

## External challenge sources

- Ota et al., *Cell* 2021, doi:10.1016/j.cell.2021.03.056, NBDC E-GEAD-397.
- Höllbacher et al., *ImmunoHorizons* 2020,
  doi:10.4049/immunohorizons.2000037, GEO GSE149090.
- Norman et al., *Science* 2019, doi:10.1126/science.aax4438, GEO GSE133344.

The committed case-study target does not re-download these studies independently. It uses
the Zhu supplementary polarization table and the frozen upstream transformation in
`src/4_polarization_signatures/polarization_signature.ipynb`: same-sign Ota/Höllbacher
Wald z-scores are averaged, sign-flipped from Th2-vs-Th1 to the Th1-like direction, and
intersected with screen genes. See [`docs/METHODS.md`](../docs/METHODS.md).

## Data policy

- Never force-add H5AD, H5MU, CSV, NPZ, or raw-count payloads under `data/`.
- The current case study uses a donor-collapsed differential-expression summary, not raw
  cells or donor-level effects.
- Selected frozen result tables live in [`results/evidence/`](../results/evidence/).
- A future real-data runner must accept an explicit external input path, record its hash,
  and fail closed when the file is unavailable.
