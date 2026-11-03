from processing.atac import load_atac
import anndata as ad
import pandas as pd
from processing.labels import add_labels, clean, mark_shared
from processing.rna import load_droplet, load_facs

from conftest import make_adata


def test_clean():
    assert clean(" Heart-and Aorta_") == "heart_and_aorta_"


def test_add_labels_applies_crosswalk_to_both_atlases():
    labels = pd.DataFrame(
        {"tissue": ["Marrow", "BoneMarrow"], "cell_type": ["naive B cell", "B cells"]}, index=["tm", "cus"]
    )
    adata = ad.AnnData(obs=pd.DataFrame(index=["tm", "cus", "missing"]))
    add_labels(adata, labels)
    assert adata.obs["tissue"].tolist() == ["bonemarrow", "bonemarrow", "unknown"]
    assert adata.obs["cell_type"].tolist() == ["b_cells", "b_cells", "unknown"]


def test_load_facs_concatenates_tissues_and_drops_unannotated(facs):
    adata = load_facs(*facs)
    assert sorted(adata.obs_names) == ["rna:facs:A1.x", "rna:facs:B1.x", "rna:facs:C1.x"]
    assert adata.obs.loc["rna:facs:A1.x", ["tissue", "cell_type"]].tolist() == ["liver", "hepatocyte"]
    assert adata.var_names.tolist() == ["GENE1", "GENE2"]


def test_load_droplet_keys_cells_by_channel_and_drops_unannotated(droplet):
    adata = load_droplet(*droplet)
    assert sorted(adata.obs_names) == ["rna:droplet:10X_P1_1_AAAA", "rna:droplet:10X_P1_2_AAAA"]
    assert set(adata.obs["cell_type"]) == {"hepatocyte"}


def test_load_atac_transposes_and_labels(atac_files):
    adata = load_atac(*atac_files)
    assert adata.shape == (2, 2)
    assert adata.obs_names.tolist() == ["atac:c1", "atac:c2"]
    assert adata.obs["cell_type"].tolist() == ["hepatocyte", "unknown"]


def test_mark_shared_ignores_unknown_cell_types():
    rna = make_adata("rna", ["liver", "heart", "lung"], ["hepatocyte", "unknown", "t_cell"])
    atac = make_adata("atac", ["liver", "heart"], ["hepatocyte", "unknown"])
    mark_shared(rna, atac)
    assert rna.obs["scjoint_eligible"].tolist() == [True, False, False]
    assert atac.obs["scjoint_eligible"].tolist() == [True, False]
