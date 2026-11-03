"""Tabula Muris RNA data

These are generated from two technologies, Smart-seq2 (FACS) and 10x (droplet).
We drop cells that are missing from the annotation tables.
"""

import anndata as ad
import pandas as pd
import scanpy as sc
from scipy import sparse

from .labels import add_labels


def read_labels(csv):
    labels = pd.read_csv(csv, index_col="cell", usecols=["cell", "tissue", "cell_ontology_class"])
    return labels.rename(columns={"cell_ontology_class": "cell_type"})


def load_facs(facs_dir, labels_csv):
    """Concatenate the per-tissue `FACS/*-counts.csv` files (genes x cells)."""
    # read raw sources
    labels = read_labels(labels_csv)
    tissues = []
    for csv in sorted((facs_dir / "FACS").glob("*-counts.csv")):
        tissue = sc.read_csv(csv).T
        tissue.X = sparse.csr_matrix(tissue.X)
        tissues.append(tissue[tissue.obs_names.isin(labels.index)])

    # combine and annotate
    facs = ad.concat(tissues)
    add_labels(facs, labels)
    facs.obs["protocol"] = "facs"
    facs.obs_names = "rna:facs:" + facs.obs_names
    facs.var_names = facs.var_names.str.upper()
    return facs


def load_droplet(droplet_dir, labels_csv):
    """Concatenate the 10x libraries `droplet/<tissue>-<channel>/`.

    The keys are coded as `<channel>_<barcode>`, with the suffix `-1` dropped.
    """
    # read raw sources
    labels = read_labels(labels_csv)
    libraries = []
    for mtx in sorted((droplet_dir / "droplet").glob("*/matrix.mtx")):
        channel = mtx.parent.name.rsplit("-", 1)[1]
        library = sc.read_10x_mtx(mtx.parent)
        library.obs_names = channel + "_" + library.obs_names.str.removesuffix("-1")
        libraries.append(library[library.obs_names.isin(labels.index)])

    # combine and annotate
    droplet = ad.concat(libraries)
    add_labels(droplet, labels)
    droplet.obs["protocol"] = "droplet"
    droplet.obs_names = "rna:droplet:" + droplet.obs_names
    droplet.var_names = droplet.var_names.str.upper()
    return droplet
