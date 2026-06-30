# tests/test_core.py
import numpy as np
import scanpy as sc
import pytest
from tuf import compute_tes, compute_tds, compute_tuf


@pytest.fixture
def dummy_adata():
    """Create a minimal AnnData with neighbors + pseudotime."""
    adata = sc.datasets.pbmc3k()  # or synthetic
    sc.pp.filter_genes(adata, min_counts=1)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata, n_neighbors=15)
    sc.tl.diffmap(adata)
    sc.tl.dpt(adata)
    return adata


def test_compute_tes_adds_column(dummy_adata):
    compute_tes(dummy_adata)
    assert "traj_unc_tes" in dummy_adata.obs.columns
    assert dummy_adata.obs["traj_unc_tes"].notna().all()


def test_compute_tds_adds_column(dummy_adata):
    compute_tds(dummy_adata)
    assert "traj_unc_tds" in dummy_adata.obs.columns
    vals = dummy_adata.obs["traj_unc_tds"]
    assert (vals >= 0).all() and (vals <= 1).all()


def test_compute_tuf_adds_both(dummy_adata):
    compute_tuf(dummy_adata)
    assert "traj_unc_tes" in dummy_adata.obs.columns
    assert "traj_unc_tds" in dummy_adata.obs.columns


def test_missing_pseudotime_raises():
    adata = sc.AnnData(np.random.randn(50, 20))
    with pytest.raises(KeyError):
        compute_tes(adata)


def test_missing_neighbors_raises():
    adata = sc.AnnData(np.random.randn(50, 20))
    adata.obs["pseudotime"] = np.random.rand(50)
    with pytest.raises(KeyError):
        compute_tes(adata)
