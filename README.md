# TUF: Trajectory Uncertainty Framework

[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

**TUF** is a lightweight Python package for **quantifying uncertainty in single-cell trajectory inference**. Rather than inferring trajectories itself, TUF operates downstream of existing trajectory inference methods to identify regions of developmental trajectories exhibiting **high temporal uncertainty** and **directional ambiguity**.

The framework introduces two complementary uncertainty metrics:

- **Temporal Entropy Score (TES):** quantifies local temporal mixing among neighboring cells.
- **Trajectory Divergence Score (TDS):** measures directional ambiguity and branch divergence within the local neighborhood.

TUF integrates seamlessly with **Scanpy** workflows and is compatible with any trajectory inference method that stores pseudotime values in an `AnnData` object.

---

## Why TUF?

Trajectory inference methods estimate developmental progression but generally do not quantify the confidence of those estimates. TUF complements existing trajectory inference algorithms by providing local uncertainty measures that help identify:

- Transitional cell states
- Regions with mixed temporal identity
- Branch points with ambiguous developmental direction
- Cells where trajectory assignments are less reliable

TUF is designed as a **post-hoc analysis framework**, allowing uncertainty quantification without modifying existing trajectory inference pipelines.

---

## Features

- Native support for `AnnData`
- Seamless integration with Scanpy
- Automatic pseudotime detection
- Graph-based uncertainty quantification
- Publication-quality visualization functions
- Lightweight with minimal dependencies
- Compute TES and TDS individually or together

---

## Installation

Currently, TUF can be installed directly from the GitHub repository.

```bash
git clone https://github.com/Sizerta/tuf_python.git
cd tuf_python
pip install -e .
```

PyPI support is planned for a future release.

---

## Quick Start

```python
import scanpy as sc
from tuf import compute_tuf, plot_tuf

# Load example dataset
adata = sc.datasets.pbmc3k()

# Standard Scanpy preprocessing
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)
adata = adata[:, adata.var.highly_variable]

sc.pp.pca(adata)
sc.pp.neighbors(adata)

# Compute pseudotime
sc.tl.diffmap(adata)
sc.tl.dpt(adata)

# Compute trajectory uncertainty
compute_tuf(adata)

# Visualize results
plot_tuf(adata)
```

---

## Computing Individual Metrics

Compute each uncertainty metric independently if desired.

```python
from tuf import compute_tes, compute_tds

compute_tes(adata)
compute_tds(adata)
```

Visualize the results

```python
from tuf import plot_tes, plot_tds

plot_tes(adata)
plot_tds(adata)
```

---

## Output

After computation, TUF stores the uncertainty scores in `adata.obs`.

| Column | Description |
|---------|-------------|
| `traj_unc_tes` | Temporal Entropy Score (TES) |
| `traj_unc_tds` | Trajectory Divergence Score (TDS) |

These scores can be used for downstream visualization, statistical analysis, and identification of uncertain or transitional cell populations.

---

## API

| Function | Description |
|----------|-------------|
| `compute_tes()` | Compute the Temporal Entropy Score |
| `compute_tds()` | Compute the Trajectory Divergence Score |
| `compute_tuf()` | Compute both TES and TDS |
| `plot_tes()` | Visualize TES on a UMAP embedding |
| `plot_tds()` | Visualize TDS on a UMAP embedding |
| `plot_tuf()` | Display TES, TDS, and pseudotime together |

---

## Requirements

- Python ≥ 3.10
- NumPy
- Scanpy
- AnnData
- Matplotlib

---

## Citation

If you use **TUF** in your research, please cite the GitHub repository until the accompanying manuscript is published.

```text
Masoud Mahdavifar.
TUF: Trajectory Uncertainty Framework.
GitHub repository.
https://github.com/Sizerta/tuf_python
```

The manuscript describing the methodology is currently in preparation.

---

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

## Contributing

Contributions are welcome. If you discover a bug, have a feature request, or would like to contribute improvements, please open an issue or submit a pull request.

---

## Acknowledgments

TUF is built upon the excellent scientific Python ecosystem, particularly:

- Scanpy
- AnnData
- NumPy
- Matplotlib