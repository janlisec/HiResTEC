
<!-- README.md is generated from README.Rmd. Please edit that file -->

# HiResTEC

<!-- badges: start -->

[![R-CMD-check](https://github.com/janlisec/HiResTEC/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/janlisec/HiResTEC/actions/workflows/R-CMD-check.yaml)
[![CRAN
status](https://www.r-pkg.org/badges/version/HiResTEC)](https://CRAN.R-project.org/package=HiResTEC)
[![Static
Badge](https://img.shields.io/badge/doi-10.1021/acs.analchem.8b00356-yellow.svg)](https://doi.org/10.1021/acs.analchem.8b00356)
<!-- badges: end -->

Identifying labeled compounds in a <sup>13</sup>C-tracer experiment in
non-targeted fashion is a cumbersome process. This package facilitates
such type of analyses by providing high level quality control plots,
spectra deconvolution and evaluation and performing a multitude of tests
in an automatic fashion. The main idea is to use changing intensity
ratios of ion pairs from peak list, generated with i.e. `xcms`, as
candidates and evaluate those against base peak chromatograms and
spectra information within the raw measurement data automatically. The
functionality is described in detail in the publication by Hoffmann et
al. (2018) <doi:10.1021/acs.analchem.8b00356>.

## Installation

You can install the development version of HiResTEC from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("janlisec/HiResTEC")
```
