"""Combine the RNA and ATAC modalities into one MuData."""

import muon as mu

from .scjoint import COMMIT


def build_mdata(rna, atac, checksums):
    mdata = mu.MuData({"rna": rna, "atac": atac})
    mdata.uns["provenance"] = checksums
    mdata.uns["source_studies"] = {
        "rna": "Tabula Muris FACS and droplet data for gene expression",
        "atac": "Cusanovich mouse atlas, binarized gene activity from gene accessibility",
    }
    mdata.uns["modelling_distributions"] = {"rna": "negative_binomial", "atac": "bernoulli"}
    mdata.uns["scjoint_commit"] = COMMIT
    return mdata
