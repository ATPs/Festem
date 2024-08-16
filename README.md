# Festem

<!-- badges: start -->

[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![Last
Commit](https://badgen.net/github/last-commit/Edward-Z-Chen/Festem/main)
![Commits Since
Latest](https://img.shields.io/github/commits-since/Edward-Z-Chen/Festem/latest/main)
[![R-CMD-check](https://github.com/Edward-Z-Chen/Festem/actions/workflows/main_check.yaml/badge.svg?branch=main)](https://github.com/Edward-Z-Chen/Festem/actions/workflows/main_check.yaml)
<!-- badges: end -->

Festem is a statistical method for the direct selection of cell-type markers for downstream clustering. Festem distinguishes marker genes with heterogeneous distribution across cells that are cluster informative. 

<img src="https://github.com/Edward-Z-Chen/Festem/blob/main/img/graphical_abstract.png?raw=true">

For analysis codes for paper "Directly selecting differentially expressed genes for single-cell clustering analyses", see https://github.com/XiDsLab/Festem_paper.

## Installation

To install Festem, start R and enter:

    if (!require("pak", quietly = TRUE))
        install.packages("pak")
    pak::pkg_install("XiDsLab/Festem")

## Usage

For detailed usage of Festem, please refer to our protocol at XXX.

## Citations

Chen, Z., Wang, C., Huang, S., Shi, Y., & Xi, R. (2024). Directly selecting cell-type marker genes for single-cell clustering analyses. Cell Reports Methods, 4(7). <https://doi.org/10.1016/j.crmeth.2024.100810>