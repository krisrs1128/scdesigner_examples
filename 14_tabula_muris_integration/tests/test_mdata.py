import muon as mu
import pandas as pd
import scipy.sparse
from processing.mdata import build_mdata
from conftest import make_adata


def test_mdata_round_trip_keeps_unpaired_modalities(tmp_path):
    rna = make_adata("rna", ["liver"], ["a"])
    atac = make_adata("atac", ["liver"], ["a"])
    mdata = build_mdata(rna, atac, pd.DataFrame({"source": ["x"]}))
    mdata.write(tmp_path / "out.h5mu")
    restored = mu.read(tmp_path / "out.h5mu")
    assert restored.mod["rna"].obs_names.tolist() == ["rna:0"]
    assert restored.mod["atac"].obs_names.tolist() == ["atac:0"]
    assert scipy.sparse.issparse(restored.mod["rna"].X)
    assert restored.uns["provenance"]["source"].tolist() == ["x"]
