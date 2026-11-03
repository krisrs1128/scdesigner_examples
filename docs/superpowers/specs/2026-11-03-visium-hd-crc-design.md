# Visium HD CRC input preparation design

Status: approved by the user, including the test-budget constraint below. No preprocessing implementation is included yet.

## Intent and agreed scope

Prepare one Visium HD tumor section as an AnnData input for later scdesigner-based ClusterDE analysis. The user selected **P2CRC at 8 micrometer resolution**, with an easy way to switch sections. Follow the notebook-and-small-helpers style of `14_tabula_muris_integration/`.

This is a small preprocessing addition. The implementation choices below were proposed defaults and accepted during design review.

## Approaches

- **One selected section, configurable in the notebook:** selected by the user. It supports the intended single-section analysis without mixing patient effects into spatial differences.
- Processing all three cancer sections by default: considered, but unnecessary for this illustration. The same helpers will support switching to P1CRC or P5CRC.
- Reprocessing reads, converting a normalized Seurat object, or repeating deconvolution: unnecessary because count matrices, spatial positions, and published annotations are available directly.

Use a notebook-level `SAMPLE = "P2CRC"` as the only setting the user changes. Helpers receive `sample` explicitly; they must not capture a separately defined default. A small accession lookup supports `P1CRC`, `P2CRC`, and `P5CRC`. No CLI, configuration file, workflow engine, or automatic multi-section loop.

## Sources and feasibility

Paper: [Oliveira et al., Nature Genetics, 2025](https://doi.org/10.1038/s41588-025-02193-3). The Visium HD series is [GSE280315](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE280315).

| Section | GEO accession | Count H5 download, bytes |
|---|---|---:|
| P1CRC | GSM8594567 | 122378546 |
| P2CRC | GSM8594568 | 297310187 |
| P5CRC | GSM8594569 | 239822835 |

Download only the selected section's count H5, tissue-position Parquet, scale-factor JSON, and annotation Parquet. For the first three files, use the GEO HTTPS prefix `https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8594nnn/{accession}/suppl/{accession}_{sample}_` with suffixes `filtered_feature_bc_matrix.h5`, `tissue_positions.parquet.gz`, and `scalefactors_json.json.gz`. Read the gzip-compressed auxiliary files without extracting an archive.

Use the inspected annotation release at `https://raw.githubusercontent.com/10XGenomics/HumanColonCancer_VisiumHD/df5801cf5e681b33b5026626989a7635fdf26b99/MetaData/{sample}_Metadata.parquet`. The commit pin keeps annotations reproducible.

P2CRC's downloaded H5 was inspected: **545,913 bins by 18,085 genes**, with **214,633,705 stored count entries**. All values are finite nonnegative integers. Every matrix barcode matches both spatial positions and annotation rows, and all are marked under tissue. The matrix contains **63 all-zero bins**. Gene IDs are unique; gene symbols are not (18,082 unique symbols). Supplementary Table 1 reports 545,943 bins; use the downloaded matrix's 545,913, not that published summary value.

These are unnormalized, probe-based UMI counts. Their preservation avoids a preprocessing obstacle to NB fitting; it does not establish that NB or any other marginal family fits well. Raw reads, full-resolution images, and nuclear segmentation are not needed.

## Files and responsibilities

New example directory: `15_visium_hd_crc/`.

- `15-visium-hd-crc.ipynb`: define `SAMPLE`, download inputs, call the loader, display QC and annotation coverage, write and reopen the selected section's `.h5ad`.
- `processing/__init__.py`: package marker.
- `processing/sources.py`: accession lookup, source URLs, selected-section downloads, and source SHA-256 calculation. Cache inputs in `data/raw/{sample}/`; download to temporary `.part` files and rename only after success. Reuse completed downloads. Unsupported sample names fail before network access.
- `processing/visium.py`: load sparse counts, join spatial and annotation tables by barcode, attach QC and provenance, and return AnnData. No sample-specific variable names or branches beyond the accession lookup.
- `tests/conftest.py`, `tests/test_sources.py`, and `tests/test_visium.py`: tiny local fixtures and mocked downloads; routine tests never download a full section.
- `data/processed/{sample}_008um.h5ad`: generated output, ignored by Git along with raw inputs.

Use the existing Python 3.12 scientific environment with Scanpy/AnnData, SciPy, pandas, and h5py; declare `pyarrow` for Parquet reading. Preserve unrelated environment-file contents.

## AnnData contract

- `.X`: bins by genes, sparse CSR integer counts, unchanged from the filtered Space Ranger matrix. Do not create a duplicate counts layer or `.raw` copy.
- `.obs_names`: original unique Space Ranger barcodes. One object contains one section; `.obs["sample_id"]` identifies it.
- `.var_names`: unique source gene IDs; retain source gene symbols in `.var["gene_symbols"]`. Do not merge genes with duplicate symbols or uppercase/rename them.
- `.obs`: official `in_tissue`, `array_row`, and `array_col`; `sample_id`; QC `total_counts` and `n_genes_by_counts`; and physical coordinates `x_um`, `y_um`.
- `.obsm["spatial"]`: full-resolution image coordinates in **(x, y) = (column, row)** order, in pixels. `.obsm["spatial_um"]` holds the same coordinates multiplied by the published `microns_per_pixel`; `x_um` and `y_um` expose these values as simulation covariates. Do not negate y or use the transformed X/Y columns in the annotation table.
- Preserve the published annotation fields verbatim in `.obs`: `DeconvolutionClass`, `DeconvolutionLabel1`, `DeconvolutionLabel2`, `Periphery`, `UnsupervisedL1`, `UnsupervisedL2`, `MacrophageSubtype`, and `GobletSubcluster`. Missing values remain missing and published `Unknown` labels remain distinct from missing values.
- Add `.obs["cell_type"]` as `DeconvolutionLabel1` **only for singlet bins**; otherwise leave it missing. Retain all bins regardless of annotation certainty.
- Store the scale factors under `.uns["spatial"][sample]["scalefactors"]`. Record coordinate units and axis order explicitly in `.uns["spatial_coordinate_system"]`; no image arrays are required.
- `.uns["provenance"]`: paper DOI, GEO series/accession, sample, annotation commit, source URLs/local SHA-256 values, relevant package versions, 8 micrometer bin size, and the statement that counts were not normalized or additionally filtered.

The count matrix defines the retained barcodes. Reindex positions and annotations to that exact order; never join by row position or intersect away unmatched observations. Require unique source keys and matching coordinate/annotation rows for every matrix barcode. Missing annotation *values* are allowed; missing rows are a source mismatch and must fail clearly. Check tissue flags agree for retained barcodes.

## Constraints and exclusions

- Keep the automated suite to at most 10 collected tests, comparable to the 10 tests in `14_tabula_muris_integration/tests/`. Use small shared fixtures and related assertions; do not expand parametrization, test infrastructure, or exhaustive case matrices to work around this budget. Ask before increasing it.
- Preserve all supplied filtered-matrix bins and genes, including zero-count bins/genes. Report zero-count bins in QC; leave analysis-specific filtering to the later analysis.
- No normalization, log transformation, HVG selection, imputation, random subsampling, or conversion of the full expression matrix to dense format.
- Validate that counts are finite, nonnegative integers; coordinates are finite; source bin size is 8 micrometers; and `microns_per_pixel` is finite and positive.
- Keep doublets, rejects, unannotated bins, and their original uncertainty labels. Do not interpret bin labels as guaranteed single-cell identities.
- Published cell types and tumor/periphery labels are expression-derived context, not independent ground-truth domains. Do not create a new `domain` field or infer domains at this stage.
- No ClusterDE/scdesigner fitting, marker testing, integration, deconvolution, or histology segmentation.
- Use sample-specific cache and output paths so switching `SAMPLE` cannot reuse another section's annotations or overwrite its output.

## Acceptance checks

Within the test budget, cover sample-dependent URLs and paths, partial-download safety, unsupported samples, sparse count preservation, barcode-based joins with shuffled rows, representative invalid/duplicate/missing keys, duplicate gene symbols with distinct IDs, zero-count bins, coordinate order and units, uncertain/missing annotation values, and H5AD round-trip preservation. Combine related checks in the same test rather than adding a test for every field or validation branch.

The real-data run uses P2CRC only. It must produce and reopen `15_visium_hd_crc/data/processed/P2CRC_008um.h5ad` with shape `(545913, 18085)`, unchanged sparse counts, matched finite coordinates, and preserved annotation uncertainty. Display bin/gene totals, zero-count-bin totals, and annotation coverage. Changing only `SAMPLE` must select P1CRC or P5CRC throughout the same notebook; tests exercise this without requiring additional large downloads.

The implementer may choose local helper names and test organization within these responsibilities. No additional configuration machinery is needed.
