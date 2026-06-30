import scanpy as sc
import matplotlib.pyplot as plt


def _find_pseudotime_key(adata):
    """Helper to auto-detect the pseudotime column for plotting"""
    candidates = ["dpt_pseudotime", "pseudotime", "DC1"]
    for candidate in candidates:
        if candidate in adata.obs.columns:
            return candidate
    return "pseudotime"  # Fallback


def plot_tes(adata, **kwargs):
    sc.pl.umap(adata, color="traj_unc_tes",
               title="TES - Temporal Entropy Score",
               cmap="viridis", show=False, ** kwargs)


def plot_tds(adata, **kwargs):
    sc.pl.umap(adata, color="traj_unc_tds",
               title="TDS - Trajectory Divergence Score",
               cmap="plasma", show=False,**kwargs)


def plot_tuf(adata, figsize=(18, 6), **kwargs):
    pt_key = _find_pseudotime_key(adata)  # Auto-detect the correct key!

    fig, axes = plt.subplots(1, 3, figsize=figsize)
    sc.pl.umap(adata, color="traj_unc_tes",
               ax=axes[0], title="TES", show=False, **kwargs)
    sc.pl.umap(adata, color="traj_unc_tds",
               ax=axes[1], title="TDS", show=False, **kwargs)
    sc.pl.umap(adata, color=pt_key,
               ax=axes[2], title="Pseudotime", show=False, **kwargs)
    plt.tight_layout()
    return fig
