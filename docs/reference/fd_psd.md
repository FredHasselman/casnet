# Power Spectral Density Slope (PSD).

Estimate Alpha, Hurst Exponent and Fractal Dimension through log-log
slope.

## Usage

``` r
fd_psd(
  y,
  fs = NULL,
  removeTrend = c("no", "poly", "adaptive", "bridge")[2],
  polyOrder = 1,
  standardise = c("none", "mean.sd", "median.mad")[2],
  fitMethod = c("lowest25", "Wijnants", "Hurvich-Deo")[3],
  doPlot = FALSE,
  returnPlot = FALSE,
  returnPLAW = FALSE,
  returnInfo = FALSE,
  silent = FALSE,
  noTitle = FALSE,
  tsName = "y"
)
```

## Arguments

- y:

  A numeric vector or time series object.

- fs:

  Sample rate (default = `NULL`)

- standardise:

  standardise the series (default = `TRUE`).

- fitMethod:

  Method to decide on a frequency range for log-log fit. Can be one of:
  "lowest25","Wijnants","Hurvich-Deo" (default). See details for more
  info.

- doPlot:

  Return the log-log spectrum with linear fit (default = `TRUE`).

- returnPlot:

  Return ggplot2 object (default = `FALSE`)

- returnPLAW:

  Return the power law data (default = `FALSE`)

- returnInfo:

  Return all the data used in SDA (default = `FALSE`)

- silent:

  Run in silent-ish mode (default = `TRUE)`

- noTitle:

  Do not generate a title (only the subtitle)

- tsName:

  Name of y added as a subtitle to the plot

- detrend:

  Subtract linear trend from the series (default = `TRUE`).

## Value

A list object containing:

- A data matrix `PLAW` with columns `freq.norm`, `size` and `bulk`.

- Estimate of scaling exponent `alpha` based on a fit over the lowest
  25\\

- Estimate of the the Fractal Dimension (`FD`) using conversion
  formula's reported in Hasselman(2013).

- Information output by various functions.

## Details

Calls function
[`stats::spec.pgram()`](https://rdrr.io/r/stats/spec.pgram.html) to
estimate the scaling exponent of a timeseries based on the periodogram
frequency spectrum. After detrending and normalizing the signal (if
requested),
[`stats::spec.pgram()`](https://rdrr.io/r/stats/spec.pgram.html) is
called using a cosine taper = 0.5.

A line is fitted on the periodogram in log-log coordinates. The full
range is fitted as well as one of three fit-ranges:

- `lowest25` - The 25\\

- `Wijnants` - The 50 lowest frequencies (Wijnants et al., 2012)

- `Hurvich-Deo` - The Hurvich-Deo estimate (Hurvich & Deo, 1999)

## References

Hasselman, F. (2013). When the blind curve is finite: dimension
estimation and model inference based on empirical waveforms. Frontiers
in Physiology, 4, 75. https://doi.org/10.3389/fphys.2013.00075

Hurvich, C.M., & Deo, R.R. (1999). Plug-in Selection of the Number of
Frequencies in Regression Estimates of the Memory Parameter of a Long
Memory Time Series. *Journal of Time Series Analysis, 20(3)*, 331–341.

## See also

Other Fluctuation Analyses: [`fd_RR()`](fd_RR.md),
[`fd_allan()`](fd_allan.md), [`fd_dfa()`](fd_dfa.md),
[`fd_mfdfa()`](fd_mfdfa.md), [`fd_sda()`](fd_sda.md),
[`fd_sev()`](fd_sev.md)

## Author

Fred Hasselman
