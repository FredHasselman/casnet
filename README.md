# Disclaimer:       
**Still EXPERIMENTAL after all these years...**

[![Travis-CI Build Status](https://travis-ci.org/FredHasselman/casnet.svg?branch=master)](https://travis-ci.org/FredHasselman/casnet)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/FredHasselman/casnet?branch=master&svg=true)](https://ci.appveyor.com/project/FredHasselman/casnet)
[![CRAN status](https://www.r-pkg.org/badges/version/casnet)](https://cran.r-project.org/package=casnet)
[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)


# **casnet**: An R toolbox for studying Complex Adaptive Systems and NETworks <img src="man/figures/logo.png" align="right" alt="" width="120" />

A collection of analytic tools for studying signals recorded from complex adaptive systems or networks:

* (Cross) Recurrence Quantification Analyses (Continuous and Categorical (C)RQA, Chromatic RQA, Anisotropic RQA)
* Fluctuation Analyses (DFA, PSD slope, SDA, Multifractal DFA, Wavelet Singularity Spectrum)
* Network based time series analyses and visualisation (Recurrence Networks, spiral graphs)
* Multivariate (network based) time series analyses (Dynamic Complexity, Multiplex Recurrence Networks)


## Installing **casnet**

As soon as it is published on CRAN, install the latest stable release of casnet using `utils::install.packages(casnet, dependencies = TRUE)`.

## Development version

Either use `devtools::install_github` or `remotes::install_github`

```
library(devtools)
install_github("FredHasselman/casnet")

library(remotes)
install_github("FredHasselman/casnet")
```

### Vignette build failing?

If  building the vignettes fails on installation (using `build_vignettes = TRUE`), just omit the argument and find them at online as [Articles](https://fredhasselman.com/casnet/) or locate the vignettes in the `inst/docs/` folder of the repository.
