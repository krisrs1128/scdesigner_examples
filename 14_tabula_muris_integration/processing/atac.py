"""Read the sci-ATAC gene activity scores from the Cusanovich mouse atlas."""

import anndata as ad
import pandas as pd
from scipy.io import mmread

from .labels import add_labels


def load_atac(counts, cells, genes, labels_csv):
    """Load the genes x cells Matrix Market file into an AnnData object."""
    # build the atac-seq anndata object
    atac = ad.AnnData(
        mmread(counts).T.tocsr(),
        obs=pd.DataFrame(index=cells.read_text().split()),
        var=pd.DataFrame(index=genes.read_text().split()),
    )

    # join in the cell type labels.
    labels = pd.read_csv(labels_csv, sep="\t", index_col="cell", usecols=["cell", "tissue", "cell_label"])
    add_labels(atac, labels.rename(columns={"cell_label": "cell_type"}))
    atac.obs["protocol"] = "sci_atac"
    atac.obs_names = "atac:" + atac.obs_names
    atac.var_names = atac.var_names.str.upper()
    return atac
