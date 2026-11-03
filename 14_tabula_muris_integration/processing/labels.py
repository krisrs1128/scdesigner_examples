"""Hand-code the relationships between labels from the two atlases

The tabular Muris and Cusanovich datasets annotate tissues and cell types
slightly differently. To support shared simulation, we want to have unified
annotation. Here, we just manually annotate. Note that the Cusanovich
annotations are quite a bit more coarse, which is why we have question marks for
many of them.
"""

import re

TISSUES = {
    # Tabula Muris -> shared
    "Marrow": "BoneMarrow",
    "Large_Intestine": "LargeIntestine",
    "Heart_and_Aorta": "Heart",
    "Aorta": "Heart",  # ?
    "Brain_Myeloid": "Brain",  # ? TM has no brain regions
    "Brain_Non-Myeloid": "Brain",  # ?
    # Cusanovich -> shared
    "WholeBrain": "Brain",  # ?
    "Cerebellum": "Brain",  # ?
    "PreFrontalCortex": "Brain",  # ?
}

CELL_TYPES = {
    # Tabula Muris -> Cusanovich
    "B cell": "B cells",
    "naive B cell": "B cells",
    "immature B cell": "Immature B cells",
    "precursor B cell": "Immature B cells",  # ?
    "late pro-B cell": "Immature B cells",  # ?
    "early pro-B cell": "Immature B cells",  # ?
    "Fraction A pre-pro B cell": "Immature B cells",  # ?
    "T cell": "T cells",
    "immature T cell": "T cells",  # ?
    "DN1 thymic pro-T cell": "T cells",  # ?
    "regulatory T cell": "Regulatory T cells",
    "natural killer cell": "NK cells",
    "mature natural killer cell": "NK cells",
    "immature natural killer cell": "NK cells",  # ?
    "pre-natural killer cell": "NK cells",  # ?
    "macrophage": "Macrophages",
    "Kupffer cell": "Macrophages",
    "alveolar macrophage": "Alveolar macrophages",
    "microglial cell": "Microglia",
    "dendritic cell": "Dendritic cells",
    "Langerhans cell": "Dendritic cells",  # ?
    "professional antigen presenting cell": "Dendritic cells",  # ?
    "monocyte": "Monocytes",
    "classical monocyte": "Monocytes",
    "non-classical monocyte": "Monocytes",
    "promonocyte": "Monocytes",  # ?
    "hematopoietic precursor cell": "Hematopoietic progenitors",
    "Slamf1-negative multipotent progenitor cell": "Hematopoietic progenitors",
    "Slamf1-positive multipotent progenitor cell": "Hematopoietic progenitors",
    "granulocyte monocyte progenitor cell": "Hematopoietic progenitors",
    "common lymphoid progenitor": "Hematopoietic progenitors",
    "megakaryocyte-erythroid progenitor cell": "Hematopoietic progenitors",  # ?
    "granulocytopoietic cell": "Hematopoietic progenitors",  # ?
    "proerythroblast": "Erythroblasts",
    "erythroblast": "Erythroblasts",
    "cardiac muscle cell": "Cardiomyocytes",
    "hepatocyte": "Hepatocytes",
    "type II pneumocyte": "Type II pneumocytes",
    "oligodendrocyte": "Oligodendrocytes",
    "astrocyte": "Astrocytes",
    "Bergmann glial cell": "Astrocytes",  # ?
    "kidney proximal straight tubule epithelial cell": "Proximal tubule S3",
    "epithelial cell of proximal tubule": "Proximal tubule",
    "kidney loop of Henle ascending limb epithelial cell": "Loop of henle",
    "kidney collecting duct epithelial cell": "Collecting duct",
    "enterocyte of epithelium of large intestine": "Enterocytes",
    "epithelial cell of large intestine": "Enterocytes",  # ?
    # Tabula Muris endothelial cells are not split into I and II
    "endothelial cell": "Endothelial cells",
    "lung endothelial cell": "Endothelial cells",
    "kidney capillary endothelial cell": "Endothelial cells",
    "endothelial cell of hepatic sinusoid": "Endothelial cells",
    # Cusanovich
    "Endothelial I cells": "Endothelial cells",  # ?
    "Endothelial II cells": "Endothelial cells",  # ?
    "Endothelial I (glomerular)": "Endothelial cells",  # ?
    # Left unmapped: leukocyte, lymphocyte, myeloid cell, blood cell, granulocyte,
    # basophil, mast cell, erythrocyte, neuron, stromal and mesenchymal cells,
    # and tissues Cusanovich lacks (bladder, skin, tongue, pancreas, ...).
}


def clean(label):
    """Lowercase and join with underscores: 'Heart-and Aorta' -> 'heart_and_aorta'."""
    return re.sub(r"[\s_-]+", "_", str(label).strip().casefold())


def add_labels(adata, labels):
    """Set `obs.tissue` and `obs.cell_type`. Unlisted cells are labeled 'unknown'."""
    labels = labels.reindex(adata.obs_names).fillna("unknown")
    adata.obs["tissue"] = labels["tissue"].replace(TISSUES).map(clean).to_numpy()
    adata.obs["cell_type"] = labels["cell_type"].replace(CELL_TYPES).map(clean).to_numpy()


def strata(adata):
    return list(zip(adata.obs["tissue"], adata.obs["cell_type"]))


def mark_shared(rna, atac):
    """Label `obs.scjoint_eligible` for cells tissue/type that occurs in both modalities."""
    shared = set(strata(rna)) & set(strata(atac))
    for adata in (rna, atac):
        adata.obs["scjoint_eligible"] = [
            s in shared and s[1] != "unknown" for s in strata(adata)
        ]
