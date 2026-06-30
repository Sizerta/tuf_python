import numpy as np
from anndata import AnnData



def compute_tes(adata: AnnData,
                pseudotime_key: str = None,
                neighbors_key: str = None,
                key_added: str = "traj_unc_tes") -> None:
    """Compute Temporal Entropy Score (TES)"""

    # Auto-detect pseudotime column
    if pseudotime_key is None:
        candidates = ["dpt_pseudotime", "pseudotime", "DC1", "diffmap"]
        for candidate in candidates:
            if candidate in adata.obs.columns:
                pseudotime_key = candidate
                print(f"Auto-detected pseudotime column: {pseudotime_key}")
                break
        else:
            raise KeyError(
                "No pseudotime column found. Run sc.tl.dpt() first.")

    if pseudotime_key not in adata.obs.columns:
        raise KeyError(f"Pseudotime column '{pseudotime_key}' not found.")

    pt = adata.obs[pseudotime_key].values

    # FIX: Handle Scanpy's default neighbor graph location
    if neighbors_key is None:
        g = adata.obsp.get("connectivities")
    else:
        g = adata.obsp.get(f"{neighbors_key}_connectivities")

    if g is None:
        raise KeyError(
            "Neighbor graph not found. Run sc.pp.neighbors() first.")

    n_cells = adata.n_obs
    tes = np.zeros(n_cells)

    for i in range(n_cells):
        nbrs = g[i].nonzero()[1]
        if len(nbrs) < 3:
            continue
        nbr_pt = pt[nbrs]
        k = len(nbr_pt)
        sorted_pt = np.sort(nbr_pt)
        weighted = np.sum((2 * np.arange(1, k+1) - k - 1) * sorted_pt)
        tes[i] = weighted / (k * (k - 1) / 2)

    adata.obs[key_added] = tes



def compute_tds(adata: AnnData,
                pseudotime_key: str = None,
                reduction: str = "X_pca",
                neighbors_key: str = None,
                key_added: str = "traj_unc_tds") -> None:
    """Compute Trajectory Divergence Score (TDS)"""

    if pseudotime_key is None:
        candidates = ["dpt_pseudotime", "pseudotime", "DC1"]
        for candidate in candidates:
            if candidate in adata.obs.columns:
                pseudotime_key = candidate
                break

    if pseudotime_key not in adata.obs.columns:
        raise KeyError(f"Pseudotime column '{pseudotime_key}' not found.")

    pt = adata.obs[pseudotime_key].values
    X = adata.obsm[reduction]

    # FIX: Handle Scanpy's default neighbor graph location
    if neighbors_key is None:
        g = adata.obsp.get("connectivities")
    else:
        g = adata.obsp.get(f"{neighbors_key}_connectivities")

    if g is None:
        raise KeyError("Neighbor graph not found.")

    n_cells = adata.n_obs
    tds = np.zeros(n_cells)

    for i in range(n_cells):
        nbrs = g[i].nonzero()[1]
        if len(nbrs) < 3:
            continue

        self_pt = pt[i]
        forward_mask = pt[nbrs] > self_pt

        if np.sum(forward_mask) < 3:
            continue

        forward_idx = nbrs[forward_mask]
        center = X[i]
        vecs = X[forward_idx] - center
        norms = np.linalg.norm(vecs, axis=1)
        norms[norms < 1e-12] = 1.0
        vecs /= norms[:, np.newaxis]

        mean_dir = np.mean(vecs, axis=0)
        tds[i] = 1.0 - np.linalg.norm(mean_dir)

    adata.obs[key_added] = tds



def compute_tuf(adata: AnnData, **kwargs):
    """Compute full TUF (TES + TDS)"""
    compute_tes(adata, **kwargs)
    compute_tds(adata, **kwargs)
