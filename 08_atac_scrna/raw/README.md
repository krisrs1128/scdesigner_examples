# Raw Data — ATAC + scRNA Integration

The processing script reads three input files from this directory and writes outputs to `../data/`.

| File | Role | Source |
|------|------|--------|
| `scATAC-Healthy-Hematopoiesis-191120.rds` | **Input** — scATAC-seq `SummarizedExperiment` for healthy hematopoiesis | [Granja et al. 2019](https://www.nature.com/articles/s41587-019-0332-7) |
| `scRNA-Healthy-Hematopoiesis-191120.rds` | **Input** — scRNA-seq `SingleCellExperiment` for healthy hematopoiesis | [Granja et al. 2019](https://www.nature.com/articles/s41587-019-0332-7) |
| `scATAC-scRNA-mappings.rds` | **Input** — cell-to-cell mappings between ATAC and RNA modalities | [Granja et al. 2019](https://www.nature.com/articles/s41587-019-0332-7) |
| `MPAL-Significant-Peak2Gene-Links.tsv.gz` | Supplementary peak–gene links (not used by the current processing script) | [Granja et al. 2019](https://www.nature.com/articles/s41587-019-0332-7) |

Outputs (`atac_healthy.h5ad`, `rna_healthy.h5ad`, `atac_healthy.rds`, `rna_healthy.rds`) are written to `../data/` by `processing/process.qmd`.

Place the three `.rds` input files here before running the processing script.
