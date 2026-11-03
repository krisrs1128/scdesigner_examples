"""One tiny P2CRC section in the layout of the real downloads."""

import gzip
import io
import json

import h5py
import numpy as np
import pandas as pd
import pytest
from processing.sources import source_urls
from scipy import sparse

BARCODES = [f"b{i}" for i in range(5)]
ORDER = ["b2", "b0", "b4", "b1", "b3"]
COUNTS = np.array([[1, 0, 4], [0, 0, 0], [2, 0, 0], [0, 0, 3], [1, 0, 0]])


def read_gz_parquet(path):
    return pd.read_parquet(io.BytesIO(gzip.decompress(path.read_bytes())))


def write_gz_parquet(path, table):
    buffer = io.BytesIO()
    table.to_parquet(buffer)
    path.write_bytes(gzip.compress(buffer.getvalue()))


def read_gz_json(path):
    return json.loads(gzip.decompress(path.read_bytes()))


def write_gz_json(path, values):
    path.write_bytes(gzip.compress(json.dumps(values).encode()))


def write_counts(path):
    # A gene-by-bin CSC matrix has the same arrays as a bins-by-genes CSR matrix.
    csr = sparse.csr_matrix(COUNTS.astype("int64"))
    with h5py.File(path, "w") as f:
        m = f.create_group("matrix")
        m["data"], m["indices"], m["indptr"] = csr.data, csr.indices.astype("int64"), csr.indptr.astype("int64")
        m["shape"] = np.array([3, 5], dtype="int32")
        m["barcodes"] = np.array(BARCODES, dtype="S")
        features = m.create_group("features")
        features["id"] = np.array(["g1", "g2", "g3"], dtype="S")
        features["name"] = np.array(["DUP", "DUP", "OTHER"], dtype="S")
        features["feature_type"] = np.array(["Gene Expression"] * 3, dtype="S")
        features["genome"] = np.array(["GRCh38"] * 3, dtype="S")


def positions_table():
    i = np.arange(5)
    table = pd.DataFrame(
        {
            "barcode": BARCODES,
            "in_tissue": np.ones(5, dtype="uint8"),
            "array_row": np.zeros(5, dtype="uint32"),
            "array_col": i.astype("uint32"),
            "pxl_row_in_fullres": 10.0 * (i + 1),
            "pxl_col_in_fullres": 100.0 * (i + 1),
        }
    )
    return table.set_index("barcode").loc[ORDER].reset_index()


def annotation_table():
    table = pd.DataFrame(
        {
            "barcode": BARCODES,
            "tissue": pd.Categorical(["1"] * 5),
            "X": np.arange(5) * 7.0,
            "Y": np.arange(5) * -3.0,
            "DeconvolutionClass": ["singlet", "doublet_certain", "doublet_uncertain", "reject", None],
            "DeconvolutionLabel1": ["Tumor", "CAF", "T cell", "Tumor", None],
            "DeconvolutionLabel2": [None, "Macrophage", "CAF", None, None],
            "Periphery": ["Tumor", "50 micron", "Tissue", "Tissue", "Tissue"],
            "UnsupervisedL1": ["Tumor", "Fibroblast", "T cells", "Unknown", None],
            "UnsupervisedL2": [None] * 5,
            "MacrophageSubtype": [None] * 5,
            "GobletSubcluster": [None] * 5,
        }
    )
    return table.set_index("barcode").loc[ORDER].reset_index()


@pytest.fixture
def section_files(tmp_path):
    names = {key: url.rsplit("/", 1)[1] for key, url in source_urls("P2CRC").items()}
    paths = {key: tmp_path / "P2CRC" / name for key, name in names.items()}
    paths["counts"].parent.mkdir()
    write_counts(paths["counts"])
    write_gz_parquet(paths["positions"], positions_table())
    write_gz_json(paths["scales"], {"bin_size_um": 8, "microns_per_pixel": 0.25, "spot_diameter_fullres": 32})
    annotation_table().to_parquet(paths["annotations"])
    return paths
