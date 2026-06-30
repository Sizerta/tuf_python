import scanpy as sc
import matplotlib.pyplot as plt
from tuf.core import compute_tuf
from tuf.plotting import plot_tuf, plot_tes, plot_tds

print("🚀 Testing TUF Package...")

# ====================== 1. LOAD TEST DATA ======================
print("Loading test data...")
adata = sc.read_h5ad("endocrinogenesis_day15.5_preprocessed.h5ad")
print(f"Loaded {adata.n_obs} cells × {adata.n_vars} genes")

# ====================== 2. PREPROCESS ======================
print("Preprocessing...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.pp.pca(adata, n_comps=50)
sc.pp.neighbors(adata, n_neighbors=15)

# ====================== 3. COMPUTE PSEUDOTIME ======================
print("Computing Diffusion Pseudotime...")
sc.tl.diffmap(adata)
sc.tl.dpt(adata, n_dcs=10)

# ====================== 4. COMPUTE TUF ======================
print("\nComputing TES + TDS...")

# Call compute_tuf WITHOUT the pseudotime_key argument.
# Your core.py is smart enough to auto-detect "dpt_pseudotime"!
compute_tuf(adata, neighbors_key="neighbors")

print("\n✅ TUF computation successful!")

# ====================== 5. VISUALIZATION ======================
print("\nGenerating plots...")

fig = plot_tuf(adata, figsize=(18, 6))
plt.savefig("tuf_test_plots.png", dpi=300, bbox_inches='tight')
plt.show()

# Individual plots
plot_tes(adata)
plt.savefig("tes_umap.png", dpi=300, bbox_inches='tight')

plot_tds(adata)
plt.savefig("tds_umap.png", dpi=300, bbox_inches='tight')

print("\n✅ All plots saved!")
