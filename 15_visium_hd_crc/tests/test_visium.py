import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from conftest import read_gz_json, read_gz_parquet, write_gz_json, write_gz_parquet
from processing import sources
from processing.sources import ANNOTATION_COMMIT, source_manifest
from processing.visium import load_section
from scipy import sparse
import h5py

COUNTS = [[1, 0, 4], [0, 0, 0], [2, 0, 0], [0, 0, 3], [1, 0, 0]]


def fails(files, key, edit, match):
    """Apply `edit` to one source file, expect a clear ValueError, then restore the file."""
    original = files[key].read_bytes()
    try:
        edit(files[key])
        with pytest.raises(ValueError, match=match):
            load_section(files, "P2CRC")
    finally:
        files[key].write_bytes(original)


def edit_matrix(edit):
    def apply(path):
        with h5py.File(path, "r+") as f:
            edit(f["matrix"])

    return apply


def float_data(value):
    def edit(m):
        data = m["data"][:].astype("float64")
        data[0] = value
        del m["data"]
        m["data"] = data

    return edit_matrix(edit)


def edit_positions(edit):
    def apply(path):
        table = read_gz_parquet(path)
        edit(table)
        write_gz_parquet(path, table)

    return apply


def edit_scales(**values):
    def apply(path):
        write_gz_json(path, {**read_gz_json(path), **values})

    return apply


def edit_annotations(edit):
    def apply(path):
        table = pd.read_parquet(path)
        table = edit(table)
        table.to_parquet(path)

    return apply


def test_counts_gene_ids_and_zero_bins_are_preserved(section_files):
    a = load_section(section_files, "P2CRC")
    assert sparse.isspmatrix_csr(a.X) and np.issubdtype(a.X.dtype, np.integer)
    assert (a.X != sparse.csr_matrix(COUNTS)).nnz == 0
    assert a.obs_names.tolist() == ["b0", "b1", "b2", "b3", "b4"]
    assert a.var_names.tolist() == ["g1", "g2", "g3"]
    assert a.var["gene_symbols"].tolist() == ["DUP", "DUP", "OTHER"]
    assert a.obs["total_counts"].tolist() == [5, 0, 2, 3, 1]
    assert a.obs["n_genes_by_counts"].tolist() == [2, 0, 1, 1, 1]
    assert a.raw is None and len(a.layers) == 0


def test_shuffled_spatial_and_uncertain_annotations_are_joined_by_barcode(section_files):
    a = load_section(section_files, "P2CRC")
    np.testing.assert_allclose(a.obsm["spatial"], [[100, 10], [200, 20], [300, 30], [400, 40], [500, 50]])
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


def test_invalid_source_keys_fail(section_files):
    fails(section_files, "counts", edit_matrix(lambda m: m["barcodes"].__setitem__(1, b"b0")), "barcode")
    fails(section_files, "counts", edit_matrix(lambda m: m["features/id"].__setitem__(1, b"g1")), "gene")
    fails(section_files, "positions", edit_positions(lambda t: t.__setitem__("barcode", ["b0"] * 5)), "barcode")
    fails(section_files, "annotations", edit_annotations(lambda t: t[t.barcode != "b0"]), "annotation")


def test_invalid_counts_fail_before_integer_conversion(section_files):
    for bad in (-1.0, 0.5, np.nan):
        fails(section_files, "counts", float_data(bad), "count")


def test_invalid_spatial_inputs_fail(section_files):
    fails(section_files, "positions", edit_positions(lambda t: t.__setitem__("pxl_col_in_fullres", np.nan)), "coordinate")
    fails(section_files, "scales", edit_scales(microns_per_pixel=0), "microns_per_pixel")
    fails(section_files, "scales", edit_scales(bin_size_um=16), "bin_size_um")

    def untissue(table):
        table["tissue"] = table["tissue"].astype(str)
        table.loc[table.barcode == "b0", "tissue"] = "0"
        return table

    fails(section_files, "annotations", edit_annotations(untissue), "tissue")


def test_h5ad_round_trip_preserves_annotations_and_provenance(section_files, tmp_path):
    a = load_section(section_files, "P2CRC")
    a.write_h5ad(tmp_path / "section.h5ad")
    b = ad.read_h5ad(tmp_path / "section.h5ad")
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


def test_notebook_uses_selected_sample(section_files, tmp_path, monkeypatch):
    notebook = Path(__file__).parent.parent / "15-visium-hd-crc.ipynb"
    code = [
        "".join(cell["source"])
        for cell in json.loads(notebook.read_text())["cells"]
        if cell["cell_type"] == "code"
    ]
    calls = []

    def download(raw_dir, sample):
        calls.append((raw_dir, sample))
        return section_files

    monkeypatch.setattr(sources, "download", download)
    monkeypatch.chdir(tmp_path)
    for sample in ("P2CRC", "P1CRC"):
        namespace = {"__name__": "__main__"}
        for source in code:
            source = source.replace('SAMPLE = "P2CRC"', f'SAMPLE = "{sample}"')
            exec(compile(source, str(notebook), "exec"), namespace)
    assert calls == [(Path("data/raw"), "P2CRC"), (Path("data/raw"), "P1CRC")]
    for sample in ("P2CRC", "P1CRC"):
        a = ad.read_h5ad(tmp_path / "data" / "processed" / f"{sample}_008um.h5ad")
        assert a.shape == (5, 3)
        assert a.obs["sample_id"].eq(sample).all()
        assert a.uns["provenance"]["sample"] == sample
