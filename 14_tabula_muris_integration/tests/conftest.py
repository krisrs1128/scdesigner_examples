"""Tiny source files in the layout of the real downloads."""

import gzip

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from scipy.io import mmwrite


def tm_labels(path, cells):
    rows = [{"cell": c, "tissue": "Liver", "cell_ontology_class": "Hepatocyte"} for c in cells]
    pd.DataFrame(rows).to_csv(path, index=False)


@pytest.fixture
def facs(tmp_path):
    """Two tissues, three cells; 'unlabelled' has no annotation."""
    (tmp_path / "FACS").mkdir()
    for tissue, cells in [("Liver", ["A1.x", "B1.x"]), ("Heart", ["C1.x", "unlabelled"])]:
        counts = pd.DataFrame([[1] * len(cells), [0] * len(cells)], index=["Gene1", "Gene2"], columns=cells)
        counts.to_csv(tmp_path / "FACS" / f"{tissue}-counts.csv")
    tm_labels(tmp_path / "labels.csv", ["A1.x", "B1.x", "C1.x"])
    return tmp_path, tmp_path / "labels.csv"


@pytest.fixture
def droplet(tmp_path):
    """Two libraries sharing a barcode; 'GGGG-1' has no annotation."""
    for library, barcodes in [("Liver-10X_P1_1", ["AAAA-1", "GGGG-1"]), ("Liver-10X_P1_2", ["AAAA-1"])]:
        folder = tmp_path / "droplet" / library
        folder.mkdir(parents=True)
        (folder / "barcodes.tsv").write_text("\n".join(barcodes) + "\n")
        (folder / "genes.tsv").write_text("g1\tGene1\ng2\tGene2\n")
        mmwrite(str(folder / "matrix.mtx"), sparse.csc_matrix(np.ones((2, len(barcodes)), dtype=int)))
    tm_labels(tmp_path / "labels.csv", ["10X_P1_1_AAAA", "10X_P1_2_AAAA"])
    return tmp_path, tmp_path / "labels.csv"


@pytest.fixture
def atac_files(tmp_path):
    with gzip.open(tmp_path / "counts.mtx.gz", "wb") as f:
        mmwrite(f, sparse.csc_matrix(np.array([[1, 0], [0, 1]])))
    (tmp_path / "cells.txt").write_text("c1\nc2\n")
    (tmp_path / "genes.txt").write_text("GENE1\nGENE2\n")
    pd.DataFrame(
        {"cell": ["c1", "c2"], "tissue": ["Liver", "Heart"], "cell_label": ["Hepatocyte", "unknown"]}
    ).to_csv(tmp_path / "labels.txt", sep="\t", index=False)
    return [tmp_path / n for n in ("counts.mtx.gz", "cells.txt", "genes.txt", "labels.txt")]


def make_adata(prefix, tissue, cell_type, eligible=None, genes=("GENE1", "GENE2")):
    n = len(tissue)
    obs = pd.DataFrame({"tissue": tissue, "cell_type": cell_type}, index=[f"{prefix}:{i}" for i in range(n)])
    if eligible is not None:
        obs["scjoint_eligible"] = eligible
    X = sparse.csr_matrix(np.arange(n * len(genes)).reshape(n, len(genes)) % 3)
    return ad.AnnData(X, obs=obs, var=pd.DataFrame(index=list(genes)))
