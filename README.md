# Disclaimer: **VERY BETA**

# **casnet**: An R toolbox for studying Complex Adaptive Systems and NETworks

A collection of analytic tools for studying signals recorded from complex adaptive systems or networks:

* Recurrence Quantification Analyses (CrossRQA, Categorical RQA, Chromatic RQA, Anisotropic RQA)
* Fluctuation Analyses (DFA varieties, PSD slope, SDA, Multifractal DFA, Wavelet Singularity Spectrum)
* Coupling Analyses (Cross Conformal Mapping, Detection of Coupling Direction, CRQA)
* Network based time series analyses (Recurrence Networks, Multifractal Spectrum Networks)


## Installing **casnet**

Either use `devtools::install_github` or download the tar/zip

### Use `devtools`

```{r}
library(devtools)
install_github("FredHasselman/casnet", build_vignettes = TRUE)
```

### Use a tar/zip build

1. Go to [`FredHasselman/casnet/pkg`](https://github.com/FredHasselman/casnet/tree/master/pkg) 
2. Select the file you need
    * `.tar.gz` is the `source` package for Mac OS or Linux.
    * `.zip` is the `source` package for Windows.
    * `.tgz` is the `binary` package for any OS, select it if you want to compile the package on your machine (you'll need a compiler see [the documentation on CRAN](https://cran.r-project.org/index.html)).
3. Install the package locally either through the RStudio GUI, or by calling `install.packages()` or `devtools::install()`


