# Visium HD CRC Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Save one configurable 8 micrometer Visium HD cancer section, defaulting to P2CRC, as an annotated, count-preserving AnnData object for later ClusterDE analysis.

**Architecture:** A notebook defines `SAMPLE` once and calls small download and loading modules. The loader constructs sparse AnnData directly from the 10x H5, joins published spatial/annotation tables by barcode, and records provenance. No analysis or multi-section workflow is introduced.

**Tech Stack:** Python 3.12, h5py, NumPy, SciPy, pandas, PyArrow, AnnData, Jupyter, pytest; reuse the existing scientific environment.

**Spec:** `docs/superpowers/specs/2026-11-03-visium-hd-crc-design.md` (approved; read it alongside this plan).

## Global Constraints

- Keep the automated suite to at most 10 collected tests, comparable to the 10 tests in `14_tabula_muris_integration/tests/`. Use small shared fixtures and related assertions; do not expand parametrization, test infrastructure, or exhaustive case matrices to work around this budget. Ask before increasing it.
- Preserve all supplied filtered-matrix bins and genes, including zero-count bins/genes. Report zero-count bins in QC; leave analysis-specific filtering to the later analysis.
- No normalization, log transformation, HVG selection, imputation, random subsampling, or conversion of the full expression matrix to dense format.
- Validate that counts are finite, nonnegative integers; coordinates are finite; source bin size is 8 micrometers; and `microns_per_pixel` is finite and positive.
- Keep doublets, rejects, unannotated bins, and their original uncertainty labels. Do not interpret bin labels as guaranteed single-cell identities.
- Published cell types and tumor/periphery labels are expression-derived context, not independent ground-truth domains. Do not create a new `domain` field or infer domains at this stage.
- No ClusterDE/scdesigner fitting, marker testing, integration, deconvolution, or histology segmentation.
- Use sample-specific cache and output paths so switching `SAMPLE` cannot reuse another section's annotations or overwrite its output.
- Use a notebook-level `SAMPLE = "P2CRC"` as the only setting the user changes. Helpers receive `sample` explicitly; they must not capture a separately defined default.

## Review Focus

- Switching sections after a cached run must select different URLs, input directories, and output names; Task 1's cache test and Task 3's notebook smoke test cover this.
- Interrupted downloads must not appear complete on the next run; Task 1 tests failure followed by retry.
- Shuffled, missing, or duplicated barcode keys must never silently misalign observations; Task 2 tests reordered joins and representative key errors.
- Image row/column order and scale factors must produce documented physical coordinates; Task 2 checks exact pixel/micrometer values and invalid spatial inputs.
- Missing/uncertain labels, zero-count bins, and repeated gene symbols must survive serialization without dropping observations or merging genes; Task 2 checks preservation and H5AD round trips.

---

## File map and execution notes

All new example paths are under `15_visium_hd_crc/`:

| File | Responsibility |
|---|---|
| `15-visium-hd-crc.ipynb` | Select the sample, invoke helpers, show QC, save/reopen H5AD |
| `processing/__init__.py` | Package marker |
| `processing/sources.py` | Source lookup, cached downloads, source hashes |
| `processing/visium.py` | Sparse loading, validated barcode joins, AnnData metadata |
| `tests/conftest.py` | One tiny on-disk section fixture |
| `tests/test_sources.py` | 3 source/download tests |
| `tests/test_visium.py` | 6 loader tests plus 1 notebook smoke test |

Modify root `environment.yml` only to add `pyarrow`. It is currently **untracked user content**: preserve it, copy it into an isolated workspace if needed, and leave it untracked unless the user approves tracking the whole file. The existing `.gitignore` already excludes example data; verify this rather than rewriting ignore rules. Do not change example 14 or unrelated untracked files.

Commands below start at the repository root. Set `PY=/opt/miniconda3/envs/scdesigner/bin/python` on this machine, or select an equivalent Python 3.12 environment. That interpreter was checked for the required packages; example 14 collected **10 tests** with it. Run example-specific suites in separate processes because both examples use a top-level `processing` package. Use an isolated workspace at execution time per the worktree skill; no implementation starts until the plan is reviewed.

## Task 1: Configurable source downloads

**Files:**
- Create: `15_visium_hd_crc/processing/__init__.py`
- Create: `15_visium_hd_crc/processing/sources.py`
- Test: `15_visium_hd_crc/tests/test_sources.py`

**Interfaces:**
- Consumes: the exact GEO URL template, accession table, and pinned annotation URL in the spec.
- Produces: `SAMPLES: dict[str, str]`, mapping `P1CRC/P2CRC/P5CRC` to `GSM8594567/GSM8594568/GSM8594569`; `ANNOTATION_COMMIT = "df5801cf5e681b33b5026626989a7635fdf26b99"`.
- Produces: `source_urls(sample: str) -> dict[str, str]`, with keys `counts`, `positions`, `scales`, `annotations`.
- Produces: `download(raw_dir: Path, sample: str) -> dict[str, Path]`, saving original URL basenames under `raw_dir / sample`.
- Produces: `source_manifest(paths: dict[str, Path], sample: str) -> dict[str, dict[str, str]]`; each source entry contains `url` and the local file's `sha256` hex digest.

- [ ] **Step 1: Write the three failing source tests.**

Mock `processing.sources.urlretrieve`; no HTTP calls. The successful stub writes `b"fixture"` to its destination and records the URL. Use these test names and assertions:

```python
# test_sources_follow_sample_and_reuse_cache
# For each of the three SAMPLES, call download(tmp_path, sample) twice.
assert set(paths) == {"counts", "positions", "scales", "annotations"}
assert all(p.parent == tmp_path / sample for p in paths.values())
assert source_urls(sample)["counts"].endswith(
    f"/{SAMPLES[sample]}_{sample}_filtered_feature_bc_matrix.h5"
)
assert f"/{ANNOTATION_COMMIT}/MetaData/{sample}_Metadata.parquet" in source_urls(sample)["annotations"]
assert len(calls) == 4  # reset the call list for each sample; second download is cached
assert source_manifest(paths, sample)["counts"] == {
    "url": source_urls(sample)["counts"],
    "sha256": hashlib.sha256(b"fixture").hexdigest(),
}

# test_invalid_sample_fails_before_download
with pytest.raises(ValueError, match="sample"):
    download(tmp_path, "P9CRC")
assert calls == []
assert list(tmp_path.iterdir()) == []

# test_interrupted_download_can_be_retried
# First stub writes partial bytes, then raises OSError; catch that error.
assert not final_path.exists()  # final_path is the expected counts path
# Replace with the successful stub and retry download(tmp_path, "P2CRC").
assert paths["counts"].read_bytes() == b"fixture"
assert not list((tmp_path / "P2CRC").glob("*.part"))
```

- [ ] **Step 2: Verify the tests fail before implementation.**

Run: `(cd 15_visium_hd_crc && "$PY" -m pytest tests/test_sources.py -q)`.
Expected: import failure for the missing `processing.sources` module/API, not a dependency failure.

- [ ] **Step 3: Implement the Task 1 interfaces in `processing/sources.py`.**

Use `urllib.request.urlretrieve`, `Path`, and `hashlib.file_digest`. Validate the sample before creating directories; reuse existing final files. Download missing files to sibling `.part` paths, remove a failed partial file, and rename only after success. Hash files as streams. Keep the accession lookup as the only sample-specific branching; introduce no configuration classes or generic downloader abstraction.

- [ ] **Step 4: Verify the source tests and data ignore rules.**

Run the Step 2 command again: expected **3 passed**. Run `git check-ignore 15_visium_hd_crc/data/raw/P2CRC/example.h5 15_visium_hd_crc/data/processed/P2CRC_008um.h5ad`: expected both paths printed.

- [ ] **Step 5: Commit the source deliverable.**

```bash
git add 15_visium_hd_crc/processing/__init__.py 15_visium_hd_crc/processing/sources.py 15_visium_hd_crc/tests/test_sources.py
git commit -m "feat: add configurable Visium HD source downloads"
```

## Task 2: Count-preserving spatial AnnData loader

**Files:**
- Create: `15_visium_hd_crc/processing/visium.py`
- Create: `15_visium_hd_crc/tests/conftest.py`
- Test: `15_visium_hd_crc/tests/test_visium.py`
- Modify: `environment.yml`, dependency list next to `pandas` (add `pyarrow` only).

**Interfaces:**
- Consumes: Task 1's paths mapping, `SAMPLES`, `ANNOTATION_COMMIT`, and `source_manifest(paths, sample)`.
- Produces: `load_section(paths: dict[str, Path], sample: str) -> anndata.AnnData`, implementing the spec's entire AnnData contract. Invalid counts, source keys, and spatial inputs raise `ValueError` with a meaningful field/category in the message; missing files retain their normal file errors.

- [ ] **Step 1: Write the shared fixture and six failing loader tests.**

Fixture `section_files` returns the Task 1 paths mapping for a tiny P2CRC section. Write an actual 10x H5 with `/matrix/{data,indices,indptr,shape,barcodes}` and `/matrix/features/{id,name,feature_type,genome}`. Use barcodes `b0` through `b4`, gene IDs `g1,g2,g3`, symbols `DUP,DUP,OTHER`, and bins-by-genes counts `[[1,0,4],[0,0,0],[2,0,0],[0,0,3],[1,0,0]]`. Store the transposed matrix as gene-by-bin CSC; use `Gene Expression` feature types and `GRCh38` genome values.

Write real gzip-wrapped position Parquet and scale JSON, plus ordinary annotation Parquet. For barcode `bi`, set full-resolution column to `100*(i+1)`, row to `10*(i+1)`, array column to `i`, array row to zero, and tissue to 1. Scale values: `bin_size_um=8`, `microns_per_pixel=0.25`, `spot_diameter_fullres=32`. Both tables use row order `[b2,b0,b4,b1,b3]`. Annotation classes in barcode order are `singlet,doublet_certain,doublet_uncertain,reject,missing`; label 1 is `Tumor,CAF,T cell,Tumor,missing`, and label 2 is missing except `Macrophage` for b1 and `CAF` for b2. Set `UnsupervisedL1` to `Tumor,Fibroblast,T cells,Unknown,missing`, `Periphery` to `Tumor,50 micron,Tissue,Tissue,Tissue`, and the other published fields to missing. Store annotation `tissue` as strings, as in the release.

Use these tests and assertions; imports are conventional and omitted here:

```python
# test_counts_gene_ids_and_zero_bins_are_preserved
a = load_section(section_files, "P2CRC")
assert sparse.isspmatrix_csr(a.X) and np.issubdtype(a.X.dtype, np.integer)
assert (a.X != sparse.csr_matrix([[1,0,4],[0,0,0],[2,0,0],[0,0,3],[1,0,0]])).nnz == 0
assert a.obs_names.tolist() == ["b0", "b1", "b2", "b3", "b4"]
assert a.var_names.tolist() == ["g1", "g2", "g3"]
assert a.var["gene_symbols"].tolist() == ["DUP", "DUP", "OTHER"]
assert a.obs["total_counts"].tolist() == [5, 0, 2, 3, 1]
assert a.obs["n_genes_by_counts"].tolist() == [2, 0, 1, 1, 1]
assert a.raw is None and len(a.layers) == 0

# test_shuffled_spatial_and_uncertain_annotations_are_joined_by_barcode
a = load_section(section_files, "P2CRC")
np.testing.assert_allclose(a.obsm["spatial"], [[100,10],[200,20],[300,30],[400,40],[500,50]])
np.testing.assert_allclose(a.obsm["spatial_um"], a.obsm["spatial"] * 0.25)
np.testing.assert_allclose(a.obs[["x_um", "y_um"]], a.obsm["spatial_um"])
assert a.obs["cell_type"].iloc[0] == "Tumor"
assert a.obs["cell_type"].iloc[1:].isna().all()
assert a.obs["DeconvolutionClass"].iloc[:4].tolist() == ["singlet", "doublet_certain", "doublet_uncertain", "reject"]
assert a.obs.loc["b1", "DeconvolutionLabel2"] == "Macrophage"
assert a.obs.loc["b1", "Periphery"] == "50 micron"
assert a.obs.loc["b3", "UnsupervisedL1"] == "Unknown"
assert pd.isna(a.obs.loc["b4", "UnsupervisedL1"])
assert a.obs["sample_id"].tolist() == ["P2CRC"] * 5

# test_invalid_source_keys_fail
# Independently mutate/restore: duplicate H5 barcode; duplicate gene ID;
# duplicate positions barcode; missing annotation row b0.
with pytest.raises(ValueError, match="barcode|gene|annotation"):
    load_section(section_files, "P2CRC")

# test_invalid_counts_fail_before_integer_conversion
# Independently replace the H5 data dataset with a floating copy whose first
# entry is -1, 0.5, or NaN. Do not assign floats into the original integer dataset.
with pytest.raises(ValueError, match="count"):
    load_section(section_files, "P2CRC")

# test_invalid_spatial_inputs_fail
# Independently mutate/restore: coordinate NaN; microns_per_pixel=0;
# bin_size_um=16; annotation tissue for b0="0".
with pytest.raises(ValueError, match="coordinate|microns_per_pixel|bin_size_um|tissue"):
    load_section(section_files, "P2CRC")

# test_h5ad_round_trip_preserves_annotations_and_provenance
# Load a, write tmp_path / "section.h5ad", and read it back as b.
assert (a.X != b.X).nnz == 0
pd.testing.assert_frame_equal(a.obs, b.obs)
pd.testing.assert_frame_equal(a.var, b.var)
np.testing.assert_allclose(a.obsm["spatial_um"], b.obsm["spatial_um"])
p = b.uns["provenance"]
assert p["geo_accession"] == "GSM8594568" and p["sample"] == "P2CRC"
assert p["annotation_commit"] == ANNOTATION_COMMIT
assert p["sources"] == source_manifest(section_files, "P2CRC")
assert not p["counts_normalized"] and p["additional_filtering"] == "none"
assert set(p["package_versions"]) == {"anndata", "numpy", "scipy", "pandas", "h5py", "pyarrow"}
```

The validation tests deliberately group a few representative failures. Use small local mutations and restore the fixture between checks; do not build a case-generation system or add parametrized tests.

- [ ] **Step 2: Verify the loader tests fail before implementation.**

Run: `(cd 15_visium_hd_crc && "$PY" -m pytest tests/test_visium.py -q)`.
Expected: import failure for missing `load_section`, not missing dependencies.

- [ ] **Step 3: Implement `load_section` in `processing/visium.py` and declare PyArrow.**

Read with h5py to preserve integer precision and avoid Scanpy's float conversion. A gene-by-bin CSC matrix's `(data, indices, indptr)` arrays directly represent its transpose as bins-by-genes CSR; use the reversed source shape without a dense intermediate. Validate counts before any integer conversion, preserve integer inputs' dtype, and compute QC from sparse counts. Require nonempty unique gene IDs/barcodes and the spec's matching source rows/tissue flags. Parse gzip wrappers explicitly and use `pandas.read_parquet` on the decompressed stream. Ignore extra table rows, reindex to matrix order, and retain missing annotation values as categorical missing values suitable for H5AD serialization.

Implement the spec's exact `.obs`, `.var`, `.obsm`, and `.uns["spatial"]` names. Define `.uns["spatial_coordinate_system"]` as `{"spatial": "full-resolution image pixels", "spatial_um": "micrometers", "axis_order": ["x", "y"], "y_direction": "down"}`. Provenance keys are `paper_doi`, `geo_series`, `geo_accession`, `sample`, `annotation_commit`, `sources`, `package_versions`, `bin_size_um`, `counts_normalized`, and `additional_filtering`; values follow the spec and tests, with `paper_doi="10.1038/s41588-025-02193-3"`, `geo_series="GSE280315"`, and `bin_size_um=8`. Obtain package versions through `importlib.metadata.version`.

Add `pyarrow` to the existing environment dependency list without changing other entries. No separate schema module, writer abstraction, or generic validation library.

- [ ] **Step 4: Verify the complete suite so far.**

Run: `(cd 15_visium_hd_crc && "$PY" -m pytest tests -q)`.
Expected: **9 passed**, with six loader tests and three source tests.

- [ ] **Step 5: Commit the loader deliverable.**

```bash
git add 15_visium_hd_crc/processing/visium.py 15_visium_hd_crc/tests/conftest.py 15_visium_hd_crc/tests/test_visium.py
git commit -m "feat: preserve Visium HD counts and spatial annotations in AnnData"
```

Leave the pre-existing untracked `environment.yml` untracked; report the added dependency at handoff.

## Task 3: Single-setting notebook and P2CRC preparation

**Files:**
- Create: `15_visium_hd_crc/15-visium-hd-crc.ipynb`
- Test: `15_visium_hd_crc/tests/test_visium.py` (add the tenth and final test).
- Generate, do not commit: `15_visium_hd_crc/data/processed/P2CRC_008um.h5ad`.

**Interfaces:**
- Consumes: `download(raw_dir, sample)` and `load_section(paths, sample)` from Tasks 1–2.
- Produces: a notebook with `SAMPLE = "P2CRC"`, `raw = Path("data/raw")`, and `processed = Path("data/processed")`; output is `processed / f"{SAMPLE}_008um.h5ad"`. No additional helper API.

- [ ] **Step 1: Write `test_notebook_uses_selected_sample` as a failing offline smoke test.**

Read the notebook through `json`, capture its code-cell sources, and execute them in a temporary working directory with `exec(compile(...))`. Mock `processing.sources.download` to record `(raw_dir, sample)` and return `section_files`. Execute twice, once unchanged and once replacing only `SAMPLE = "P2CRC"` with `SAMPLE = "P1CRC"`, using a fresh namespace each time. This is one small test for the user's configuration requirement, not a notebook testing subsystem. Expected assertions:

```python
assert calls == [(Path("data/raw"), "P2CRC"), (Path("data/raw"), "P1CRC")]
for sample in ("P2CRC", "P1CRC"):
    output = tmp_path / "data" / "processed" / f"{sample}_008um.h5ad"
    a = ad.read_h5ad(output)
    assert a.shape == (5, 3)
    assert a.obs["sample_id"].eq(sample).all()
    assert a.uns["provenance"]["sample"] == sample
```

- [ ] **Step 2: Verify the smoke test fails.**

Run: `(cd 15_visium_hd_crc && "$PY" -m pytest tests/test_visium.py::test_notebook_uses_selected_sample -q)`.
Expected: **FAIL** because the notebook file does not yet exist.

- [ ] **Step 3: Create the notebook with the Task 3 interface.**

Use short markdown explanations and ordinary Python cells, without magics or IPython-only calls. Cite the paper and sources, explain bin-level uncertainty and unchanged counts, then select the sample, download/load, print shape and QC/annotation coverage, and write the H5AD with gzip compression. Release the in-memory object before reopening with `ad.read_h5ad(..., backed="r")`, report the persisted shape, and close the backing file. Do not load a second full matrix for a display check. Keep all sample-specific output paths derived from `SAMPLE`, and preserve prior sections' files when switching.

- [ ] **Step 4: Verify the test budget and notebook smoke test.**

Run `(cd 15_visium_hd_crc && "$PY" -m pytest tests --collect-only -q)` and then the same command without `--collect-only`.
Expected: **10 tests collected**, then **10 passed**. Do not add new tests merely to exercise every validation branch.

- [ ] **Step 5: Run the notebook on P2CRC.**

From the repository root, execute the same plain-Python notebook cells with the verified interpreter:

```bash
(cd 15_visium_hd_crc && "$PY" - <<'PY'
import json
from pathlib import Path
notebook = Path("15-visium-hd-crc.ipynb")
ns = {"__name__": "__main__"}
for cell in json.loads(notebook.read_text())["cells"]:
    if cell["cell_type"] == "code":
        exec(compile("".join(cell["source"]), str(notebook), "exec"), ns)
PY
)
```

Expected: `(545913, 18085)`, **63 zero-count bins**, **322327 singlet bins**, and a saved/reopened `data/processed/P2CRC_008um.h5ad`. This is the one manual full-section acceptance run, not an additional automated integration suite. Source count arrays alone occupy about 1.7 GB; keep them sparse and avoid duplicate count layers. Reuse inspected `/tmp/scdesigner-visium-P2-*` downloads only after checking they correspond to the exact URLs and original filenames; otherwise let the notebook download them normally.

- [ ] **Step 6: Verify the saved full-section counts without reloading two AnnData objects.**

The CSR-transpose construction preserves source array order. Run:

```bash
(cd 15_visium_hd_crc && "$PY" - <<'PY'
import h5py
import numpy as np
src = "data/raw/P2CRC/GSM8594568_P2CRC_filtered_feature_bc_matrix.h5"
dst = "data/processed/P2CRC_008um.h5ad"
with h5py.File(src) as source, h5py.File(dst) as saved:
    assert tuple(saved["X"].attrs["shape"]) == (545913, 18085)
    assert saved["X/data"].shape == (214633705,)
    assert np.issubdtype(saved["X/data"].dtype, np.integer)
    for key in ("data", "indices", "indptr"):
        left, right = source[f"matrix/{key}"], saved[f"X/{key}"]
        assert left.shape == right.shape
        for start in range(0, len(left), 1_000_000):
            np.testing.assert_array_equal(left[start:start+1_000_000], right[start:start+1_000_000])
print("P2CRC sparse counts preserved exactly")
PY
)
```

Expected: `P2CRC sparse counts preserved exactly`. The backed reopen in Step 5 checks that the object is readable as AnnData; fixture tests cover coordinate/annotation persistence.

- [ ] **Step 7: Commit the notebook deliverable.**

```bash
git add 15_visium_hd_crc/15-visium-hd-crc.ipynb 15_visium_hd_crc/tests/test_visium.py
git diff --cached --check
git commit -m "feat: add configurable Visium HD preparation notebook"
```

At handoff, report the H5AD path, observed dimensions, count-preservation check, and **10-test** result. Mention the local environment change, and do not claim any ClusterDE analysis or model-fit validation has been performed.
