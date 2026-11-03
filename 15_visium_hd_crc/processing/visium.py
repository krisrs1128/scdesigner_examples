"""Load a Visium HD section as an AnnData object.

We've deliberately avoided any preprocessing so that we can model these as
negative binomial.
"""

import gzip
import io
import json
from importlib.metadata import version

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from .sources import ANNOTATION_COMMIT, SAMPLES, source_manifest

ANNOTATION_FIELDS = [
    "DeconvolutionClass",
    "DeconvolutionLabel1",
    "DeconvolutionLabel2",
    "Periphery",
    "UnsupervisedL1",
    "UnsupervisedL2",
    "MacrophageSubtype",
    "GobletSubcluster",
]
PACKAGES = ["anndata", "numpy", "scipy", "pandas", "h5py", "pyarrow"]


def read_counts(path):
    """The raw counts are stored as CSC matrices."""
    with h5py.File(path) as f:
        m = f["matrix"]
        n_genes, n_bins = (int(n) for n in m["shape"][:])
        barcodes = [b.decode() for b in m["barcodes"][:]]
        gene_ids = [g.decode() for g in m["features/id"][:]]
        symbols = [g.decode() for g in m["features/name"][:]]
        data, indices, indptr = m["data"][:], m["indices"][:], m["indptr"][:]

    if data.dtype.kind == "f":
        if not np.isfinite(data).all() or (data != np.floor(data)).any():
            raise ValueError("counts must be finite integers")
    elif data.dtype.kind not in "iu":
        raise ValueError(f"counts have unsupported dtype {data.dtype}")
    if len(data) and data.min() < 0:
        raise ValueError("counts must be nonnegative")
    if data.dtype.kind == "f":
        data = data.astype("int64")
    X = sparse.csr_matrix((data, indices, indptr), shape=(n_bins, n_genes))
    return X, barcodes, gene_ids, symbols


def align(table, barcodes, name):
    """Match the cell barcodes in the source table with the matrix counts"""
    table = table.set_index("barcode")
    missing = pd.Index(barcodes).difference(table.index)
    if len(missing):
        raise ValueError(f"{name} table lacks rows for {len(missing)} matrix barcodes, e.g. {missing[0]!r}")
    return table.loc[barcodes]


def load_section(paths, sample):
    X, barcodes, gene_ids, symbols = read_counts(paths["counts"])
    X.eliminate_zeros()

    # Read resolution information.
    with gzip.open(paths["scales"]) as f:
        scales = json.load(f)
    microns_per_pixel = scales.get("microns_per_pixel")
    if scales.get("bin_size_um") != 8:
        raise ValueError(f"bin_size_um must be 8, found {scales.get('bin_size_um')!r}")

    # Read this opposition.
    with gzip.open(paths["positions"]) as f:
        positions = pd.read_parquet(io.BytesIO(f.read()))
    positions = align(positions, barcodes, "positions")
    annotations = align(pd.read_parquet(paths["annotations"]), barcodes, "annotation")
    spatial = positions[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy(dtype=float)
    spatial_um = spatial * microns_per_pixel

    # Combined into an obs dataset
    obs = pd.DataFrame(
        {
            "in_tissue": positions["in_tissue"].to_numpy(),
            "array_row": positions["array_row"].to_numpy(),
            "array_col": positions["array_col"].to_numpy(),
            "sample_id": pd.Categorical([sample] * len(barcodes)),
            "total_counts": np.asarray(X.sum(axis=1)).ravel(),
            "n_genes_by_counts": X.getnnz(axis=1),
            "x_um": spatial_um[:, 0],
            "y_um": spatial_um[:, 1],
        },
        index=pd.Index(barcodes),
    )
    for field in ANNOTATION_FIELDS:
        obs[field] = annotations[field].astype("category")
    singlets = annotations["DeconvolutionClass"] == "singlet"
    obs["cell_type"] = annotations["DeconvolutionLabel1"].where(singlets).astype("category")

    # Build the anndata object and provide useful metadata.
    a = ad.AnnData(X=X, obs=obs, var=pd.DataFrame({"gene_symbols": symbols}, index=pd.Index(gene_ids)))
    a.obsm["spatial"] = spatial
    a.obsm["spatial_um"] = spatial_um
    a.uns["spatial"] = {sample: {"scalefactors": scales}}
    a.uns["spatial_coordinate_system"] = {
        "spatial": "full-resolution image pixels",
        "spatial_um": "micrometers",
        "axis_order": ["x", "y"],
        "y_direction": "down",
    }
    a.uns["provenance"] = {
        "paper_doi": "10.1038/s41588-025-02193-3",
        "geo_series": "GSE280315",
        "geo_accession": SAMPLES[sample],
        "sample": sample,
        "annotation_commit": ANNOTATION_COMMIT,
        "sources": source_manifest(paths, sample),
        "package_versions": {name: version(name) for name in PACKAGES},
        "bin_size_um": 8,
        "counts_normalized": False,
        "additional_filtering": "none",
    }
    return a
