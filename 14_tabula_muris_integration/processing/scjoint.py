"""Run scJoint on cells from strata shared by RNA and ATAC.

scJoint is a script, not a package: we clone a pinned commit, overwrite its
`config.py`, patch its progress bar, and run `main.py`.
"""

import subprocess
import sys

import numpy as np
import pandas as pd
from scipy import sparse

REPO_URL = "https://github.com/SydneyBioX/scJoint.git"
COMMIT = "cbbfa5d7e7faa35f3d725bd5257507a01a22810f"
N_DIMS = 64
STTY_CALL = "_, term_width = os.popen('stty size', 'r').read().split()"
TERM_WIDTH_CALL = "term_width = __import__('shutil').get_terminal_size((80, 24)).columns"

CONFIG = """\
import numpy as np
import torch

np.float = float  # scJoint still uses np.float

class Config:
    def __init__(self):
        self.use_cuda = False
        self.threads = 1
        self.device = torch.device('cpu')
        self.number_of_class = {n_types}
        self.input_size = {n_genes}
        self.rna_paths = ['{in_dir}/rna_matrix.npz']
        self.rna_labels = ['{in_dir}/rna_labels.txt']
        self.atac_paths = ['{in_dir}/atac_matrix.npz']
        self.atac_labels = []
        self.rna_protein_paths = []
        self.atac_protein_paths = []
        self.batch_size = 256
        self.lr_stage1 = 0.01
        self.lr_stage3 = 0.01
        self.lr_decay_epoch = 20
        self.epochs_stage1 = 10
        self.epochs_stage3 = 10
        self.p = 0.8
        self.embedding_size = {n_dims}
        self.momentum = 0.9
        self.center_weight = 1
        self.with_crossentorpy = True
        self.seed = 1
        self.checkpoint = ''
"""


def write_inputs(rna, atac, in_dir):
    """Save binary + integer data to run scJoint on"""
    # drop cells that aren't annotated
    in_dir.mkdir(parents=True, exist_ok=True)
    rna = rna[rna.obs["scjoint_eligible"].to_numpy()]
    atac = atac[atac.obs["scjoint_eligible"].to_numpy()]
    genes = sorted(set(rna.var_names) & set(atac.var_names))

    # save to disk
    for name, adata in [("rna", rna), ("atac", atac)]:
        counts = sparse.csr_matrix(adata[:, genes].X)
        sparse.save_npz(in_dir / f"{name}_matrix.npz", (counts > 0).astype(np.int8))
    codes, types = pd.factorize(rna.obs["cell_type"], sort=True)
    np.savetxt(in_dir / "rna_labels.txt", codes, fmt="%d")
    return len(genes), len(types)


def clone(repo):
    """Check out the pinned scJoint commit into `repo`."""
    if not repo.exists():
        subprocess.run(["git", "clone", REPO_URL, repo], check=True)
    subprocess.run(["git", "-C", repo, "checkout", "--force", COMMIT], check=True)

    # The progress bar needs to be disabled to work in a notebook
    utils = repo / "util" / "utils.py"
    utils.write_text(utils.read_text().replace(STTY_CALL, TERM_WIDTH_CALL))


def run(repo, in_dir, n_genes, n_types):
    """Train scJoint in `repo`. Returns (rna, atac) embeddings."""
    config = CONFIG.format(in_dir=in_dir.resolve(), n_genes=n_genes, n_types=n_types, n_dims=N_DIMS)
    (repo / "config.py").write_text(config)
    subprocess.run([sys.executable, "main.py"], cwd=repo, check=True)
    rna = np.loadtxt(repo / "output" / "rna_matrix_embeddings.txt", ndmin=2)
    atac = np.loadtxt(repo / "output" / "atac_matrix_embeddings.txt", ndmin=2)
    return rna, atac


def add_embeddings(adata, embedding):
    """Store the embedding in `obsm['X_scjoint']` and `obs['scjoint_1'...]`"""
    coords = np.full((adata.n_obs, N_DIMS), np.nan)
    coords[adata.obs["scjoint_eligible"].to_numpy()] = embedding
    adata.obsm["X_scjoint"] = coords
    adata.obs[[f"scjoint_{i}" for i in range(1, N_DIMS + 1)]] = coords
