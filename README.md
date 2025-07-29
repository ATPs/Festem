# FestemEnhanced

<!-- badges: start -->

[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![Last
Commit](https://badgen.net/github/last-commit/ATPs/Festem/main)
![Commits Since
Latest](https://img.shields.io/github/commits-since/ATPs/Festem/latest/main)
![GitHub Downloads](https://img.shields.io/github/downloads/ATPs/Festem/total)
<!-- badges: end -->

FestemEnhanced is an enhanced version of the [Festem package](https://github.com/XiDsLab/Festem) that adds support for normalized data analysis. This package extends the original Festem functionality by directly selecting differentially expressed genes for single-cell clustering analyses via EM-test, with added capability to work with normalized data in Seurat objects.

<img src="https://github.com/Edward-Z-Chen/Festem/blob/main/img/graphical_abstract.png?raw=true">

For analysis codes for paper "Directly selecting differentially expressed genes for single-cell clustering analyses", see https://github.com/XiDsLab/Festem_paper.

## Installation

To install FestemEnhanced, start R and enter:

    if (!require("pak", quietly = TRUE))
        install.packages("pak")
    pak::pkg_install("ATPs/Festem")

## Usage

For detailed usage of the original Festem functionality, please refer to our protocol at [STAR Protocol](https://doi.org/10.1016/j.xpro.2024.103514).

### Working with Normalized Data

For Seurat objects that only have normalized data (in the "data" layer) without raw counts, you can use the `RunFestemData` function:

```r
# Run Festem on normalized data
seurat_obj <- RunFestemData(seurat_obj)
```

This function adapts the Festem algorithm to work with normalized data by:
1. Using a modified EM test suitable for continuous normalized data
2. Skipping count-specific preprocessing steps
3. Adjusting filtering criteria for normalized values

## Key Features

- All original Festem functionality for raw count data
- New support for normalized data analysis via `RunFestemData`
- EM-test for mixture model component testing
- Integration with Seurat objects
- Parallel processing support

## Citations

- Chen, Z., Wang, C., Huang, S., Shi, Y., & Xi, R. (2024). Directly selecting cell-type marker genes for single-cell clustering analyses. Cell Reports Methods, 4(7). <https://doi.org/10.1016/j.crmeth.2024.100810>
- Wang, C., Chen, Z., & Xi, R. (2023). Feature screening for clustering analysis. arXiv preprint arXiv:2306.12671. <https://arxiv.org/abs/2306.12671>
