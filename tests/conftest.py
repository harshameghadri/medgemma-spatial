"""Shared fixtures: a small synthetic Visium-like tissue with known structure."""

import numpy as np
import pytest
import scipy.sparse as sp

T_CELL_GENES = ["CD3D", "CD3E", "TRAC", "CD2", "CD7", "PTPRC"]
EPITHELIAL_GENES = ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"]
STROMAL_GENES = ["COL1A1", "COL1A2", "FAP", "PDPN", "POSTN"]
N_BACKGROUND_GENES = 3200
GRID_SIDE = 20
SPOT_SPACING = 100.0


def make_synthetic_tissue(seed: int = 0):
    """Left half T-cell rich, right half epithelial; raw integer counts."""
    import anndata as ad

    rng = np.random.default_rng(seed)
    genes = (
        T_CELL_GENES
        + EPITHELIAL_GENES
        + STROMAL_GENES
        + ["MT-CO1"]
        + [f"BG{i}" for i in range(N_BACKGROUND_GENES)]
    )
    coords = np.array(
        [(i, j) for i in range(GRID_SIDE) for j in range(GRID_SIDE)], dtype=float
    ) * SPOT_SPACING
    n_spots = coords.shape[0]
    counts = rng.poisson(1.0, (n_spots, len(genes))).astype(np.float32)
    left = coords[:, 0] < (GRID_SIDE / 2) * SPOT_SPACING
    t_idx = [genes.index(g) for g in T_CELL_GENES]
    e_idx = [genes.index(g) for g in EPITHELIAL_GENES]
    counts[np.ix_(left, t_idx)] += rng.poisson(8, (left.sum(), len(t_idx)))
    counts[np.ix_(~left, e_idx)] += rng.poisson(8, ((~left).sum(), len(e_idx)))

    adata = ad.AnnData(sp.csr_matrix(counts))
    adata.var_names = genes
    adata.obs_names = [f"spot{i}" for i in range(n_spots)]
    adata.obsm["spatial"] = coords
    adata.obs["is_t_region"] = left
    return adata


@pytest.fixture
def synthetic_tissue():
    """Fresh synthetic tissue per test."""
    return make_synthetic_tissue()
