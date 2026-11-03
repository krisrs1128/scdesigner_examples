"""Download source files, unpack archives, record checksums."""

from hashlib import file_digest
from pathlib import Path
from urllib.request import urlretrieve
from zipfile import ZipFile

import pandas as pd

FIGSHARE = "https://ndownloader.figshare.com/files/"
CUSANOVICH = "https://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/"

# name: (local file name, url)
SOURCES = {
    "facs_counts": ("FACS.zip", FIGSHARE + "10700143"),
    "facs_labels": ("annotations_facs.csv", FIGSHARE + "13088129"),
    "droplet_counts": ("droplet.zip", FIGSHARE + "10700167"),
    "droplet_labels": ("annotations_droplet.csv", FIGSHARE + "13088039"),
    "atac_counts": (
        "activity_scores.binarized.mtx.gz",
        CUSANOVICH + "activity_score_matrices/activity_scores.binarized.mtx.gz",
    ),
    "atac_cells": (
        "activity_scores.binarized.cells.txt",
        CUSANOVICH + "activity_score_matrices/activity_scores.binarized.cells.txt",
    ),
    "atac_genes": (
        "activity_scores.binarized.genes.txt",
        CUSANOVICH + "activity_score_matrices/activity_scores.binarized.genes.txt",
    ),
    "atac_labels": ("cell_metadata.txt", CUSANOVICH + "metadata/cell_metadata.txt"),
}


def download(raw_dir):
    """Fetch missing sources into `raw_dir`; return {name: path}."""
    raw_dir = Path(raw_dir)
    raw_dir.mkdir(parents=True, exist_ok=True)
    paths = {}
    for name, (filename, url) in SOURCES.items():
        paths[name] = raw_dir / filename
        if not paths[name].exists():
            part = paths[name].with_name(filename + ".part")
            urlretrieve(url, part)
            part.rename(paths[name])
    return paths


def unzip(zip_path):
    """Extract `x.zip` to a sibling directory `x/` and return it."""
    out_dir = zip_path.with_suffix("")
    if not out_dir.exists():
        ZipFile(zip_path).extractall(out_dir)
    return out_dir