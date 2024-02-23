
<!-- README.md is generated from README.Rmd. Please edit that file -->

# HiResTEC

<!-- badges: start -->

[![R-CMD-check](https://github.com/janlisec/HiResTEC/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/janlisec/HiResTEC/actions/workflows/R-CMD-check.yaml)
[![CRAN
status](https://www.r-pkg.org/badges/version/HiResTEC)](https://CRAN.R-project.org/package=HiResTEC)
<!-- badges: end -->

Identifying labeled compounds in a <sup>13</sup>C-tracer experiment in
non-targeted fashion is a cumbersome process. This package facilitates
such type of analyses by providing high level quality control plots,
deconvoluting and evaluating spectra and performing a multitude of tests
in an automatic fashion. The main idea is to use changing intensity
ratios of ion pairs from peak list generated with ‘xcms’ as candidates
and evaluate those against base peak chromatograms and spectra
information within the raw measurement data automatically. The
functionality is described in detail in the publication by Hoffmann et
al. (2018) <doi:10.1021/acs.analchem.8b00356>.

## Installation

You can install the development version of HiResTEC from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("janlisec/HiResTEC")
```

## Example

I am sorry that I have to refer you to the publication to learn more
about HiResTEC as providing example data in the package would require to
include data in xcmsRaw format which would not be allowed on CRAN.
