"""Read the Visium HD cancer data from the links provided in

https://www.nature.com/articles/s41588-025-02193-3
"""

import hashlib
from pathlib import Path
from urllib.request import urlretrieve

SAMPLES = {"P1CRC": "GSM8594567", "P2CRC": "GSM8594568", "P5CRC": "GSM8594569"}
ANNOTATION_COMMIT = "df5801cf5e681b33b5026626989a7635fdf26b99"
GEO = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8594nnn"
ANNOTATIONS = "https://raw.githubusercontent.com/10XGenomics/HumanColonCancer_VisiumHD"


def source_urls(sample):
    if sample not in SAMPLES:
        raise ValueError(f"Unsupported sample {sample!r}; choose one of {sorted(SAMPLES)}")
    prefix = f"{GEO}/{SAMPLES[sample]}/suppl/{SAMPLES[sample]}_{sample}_"
    return {
        "counts": prefix + "filtered_feature_bc_matrix.h5",
        "positions": prefix + "tissue_positions.parquet.gz",
        "scales": prefix + "scalefactors_json.json.gz",
        "annotations": f"{ANNOTATIONS}/{ANNOTATION_COMMIT}/MetaData/{sample}_Metadata.parquet",
    }


def download(raw_dir, sample):
    urls = source_urls(sample)
    directory = Path(raw_dir) / sample
    directory.mkdir(parents=True, exist_ok=True)
    paths = {}

    # Loop over datasets
    for key, url in urls.items():
        paths[key] = final = directory / url.rsplit("/", 1)[1]
        if final.exists():
            continue
        partial = final.with_name(final.name + ".part")

        # Try downloading the data.
        try:
            urlretrieve(url, partial)
        except BaseException:
            partial.unlink(missing_ok=True)
            raise
        partial.rename(final)
    return paths


def source_manifest(paths, sample):
    urls = source_urls(sample)
    manifest = {}
    for key, path in paths.items():
        with open(path, "rb") as f:
            manifest[key] = {"url": urls[key], "sha256": hashlib.file_digest(f, "sha256").hexdigest()}
    return manifest
