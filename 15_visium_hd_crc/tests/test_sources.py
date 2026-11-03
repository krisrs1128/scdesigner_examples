import hashlib

import pytest
from processing import sources
from processing.sources import ANNOTATION_COMMIT, SAMPLES, download, source_manifest, source_urls


def succeed(calls):
    def urlretrieve(url, filename):
        calls.append(url)
        with open(filename, "wb") as f:
            f.write(b"fixture")

    return urlretrieve


def test_sources_follow_sample_and_reuse_cache(tmp_path, monkeypatch):
    for sample in SAMPLES:
        calls = []
        monkeypatch.setattr(sources, "urlretrieve", succeed(calls))
        paths = download(tmp_path, sample)
        download(tmp_path, sample)
        assert set(paths) == {"counts", "positions", "scales", "annotations"}
        assert all(p.parent == tmp_path / sample for p in paths.values())
        assert source_urls(sample)["counts"].endswith(
            f"/{SAMPLES[sample]}_{sample}_filtered_feature_bc_matrix.h5"
        )
        assert f"/{ANNOTATION_COMMIT}/MetaData/{sample}_Metadata.parquet" in source_urls(sample)["annotations"]
        assert len(calls) == 4
        assert source_manifest(paths, sample)["counts"] == {
            "url": source_urls(sample)["counts"],
            "sha256": hashlib.sha256(b"fixture").hexdigest(),
        }


def test_invalid_sample_fails_before_download(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(sources, "urlretrieve", succeed(calls))
    with pytest.raises(ValueError, match="sample"):
        download(tmp_path, "P9CRC")
    assert calls == []
    assert list(tmp_path.iterdir()) == []


def test_interrupted_download_can_be_retried(tmp_path, monkeypatch):
    def interrupted(url, filename):
        with open(filename, "wb") as f:
            f.write(b"par")
        raise OSError("connection reset")

    monkeypatch.setattr(sources, "urlretrieve", interrupted)
    final_path = tmp_path / "P2CRC" / "GSM8594568_P2CRC_filtered_feature_bc_matrix.h5"
    with pytest.raises(OSError):
        download(tmp_path, "P2CRC")
    assert not final_path.exists()

    monkeypatch.setattr(sources, "urlretrieve", succeed([]))
    paths = download(tmp_path, "P2CRC")
    assert paths["counts"].read_bytes() == b"fixture"
    assert not list((tmp_path / "P2CRC").glob("*.part"))
