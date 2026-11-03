import numpy as np
import scipy.sparse
from processing import scjoint

from conftest import make_adata


def modalities():
    rna = make_adata("rna", ["liver"] * 3, ["a", "b", "c"], eligible=[True, False, True], genes=["GENE_B", "GENE_A", "ONLY_RNA"])
    atac = make_adata("atac", ["liver"] * 2, ["a", "c"], eligible=[True, True], genes=["GENE_A", "GENE_B"])
    return rna, atac


def test_write_inputs_keeps_eligible_cells_and_shared_genes(tmp_path):
    rna, atac = modalities()
    n_genes, n_types = scjoint.write_inputs(rna, atac, tmp_path)
    assert (n_genes, n_types) == (2, 2)
    x = scipy.sparse.load_npz(tmp_path / "rna_matrix.npz")
    assert x.shape == (2, 2) and set(x.data) == {1}
    assert np.loadtxt(tmp_path / "rna_labels.txt").tolist() == [0, 1]
    assert scipy.sparse.load_npz(tmp_path / "atac_matrix.npz").shape == (2, 2)


def test_run_writes_config_and_reads_embeddings(tmp_path, monkeypatch):
    (tmp_path / "output").mkdir()
    for name in ("rna", "atac"):
        np.savetxt(tmp_path / "output" / f"{name}_matrix_embeddings.txt", np.ones((2, 64)))
    monkeypatch.setattr(scjoint.subprocess, "run", lambda *a, **k: None)
    rna, atac = scjoint.run(tmp_path, tmp_path / "inputs", n_genes=7, n_types=3)
    config = (tmp_path / "config.py").read_text()
    assert "self.input_size = 7" in config and "self.number_of_class = 3" in config
    assert rna.shape == atac.shape == (2, 64)


def test_add_embeddings_leaves_ineligible_cells_missing():
    rna, _ = modalities()
    scjoint.add_embeddings(rna, np.arange(128, dtype=float).reshape(2, 64))
    assert rna.obsm["X_scjoint"].shape == (3, 64)
    assert np.isnan(rna.obs.loc["rna:1", "scjoint_1"])
    assert rna.obs.loc["rna:2", "scjoint_2"] == 65
