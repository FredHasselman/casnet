# Standardised Dispersion Analysis (SDA).

fd_sda

## Usage

``` r
fd_sda(
  y,
  fs = NULL,
  removeTrend = c("no", "poly", "adaptive", "bridge")[2],
  polyOrder = 1,
  standardise = c("none", "mean.sd", "median.mad")[2],
  adjustSumOrder = FALSE,
  scaleMin = 4,
  scaleMax = stats::nextn(floor(NROW(y)/2), factors = 2),
  scaleResolution = log2(scaleMax) - log2(scaleMin),
  scaleS = NA,
  overlap = 0,
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

  Sample rate

- removeTrend:

  Method to use for global detrending (default = `"poly"`)

- polyOrder:

  Order of global polynomial trend to remove if `removeTrend = "poly"`.
  If `removeTrend = "adaptive"` polynomials `1` to `polyOrder` will be
  evaluated and the best fitting curve (R squared) will be removed
  (default = `1`)

- standardise:

  Standardise the series using [`ts_standardise()`](ts_standardise.md)
  with `adjustN = FALSE` (default = "mean.sd")

- adjustSumOrder:

  Adjust the time series (summation or difference), based on the global
  scaling exponent, see e.g. [Ihlen
  (2012)](https://www.frontiersin.org/files/Articles/23948/fphys-03-00141-r2/image_m/fphys-03-00141-t001.jpg)
  (default = `FALSE`)

- scaleMin:

  Minimum scale (in data points) to use for log-log regression (default
  = `4`)

- scaleMax:

  Maximum scale (in data points) to use for log-log regression (default
  = `stats::nextn(floor(NROW(y)/4), factors = 2)`)

- scaleResolution:

  The scales at which detrended fluctuation will be evaluated are
  calculated as: `seq(scaleMin, scaleMax, length.out = scaleResolution)`
  (default = `round(log2(scaleMax-scaleMin))`).

- scaleS:

  If not `NA`, it should be a numeric vector listing the scales on which
  to evaluate the detrended fluctuations. Arguments
  `scaleMax, scaleMin, scaleResolution` will be ignored (default = `NA`)

- overlap:

  A number in `[0 ... 1]` representing the amount of 'bin overlap' when
  calculating the fluctuation. This reduces impact of arbitrary time
  series begin and end points. If `length(y) = 1024` and overlap is
  `.5`, a scale of `4` will be considered a sliding window of size `4`
  with step-size `floor(.5 * 4) = 2`, so for scale `128` step-size will
  be `64` (default = `NA`)

- doPlot:

  Output the log-log scale versus fluctuation plot with linear fit by
  calling function [`plotFD_loglog()`](plotFD_loglog.md) (default =
  `TRUE`)

- returnPlot:

  Return ggplot2 object (default = `FALSE`)

- returnPLAW:

  Return the power law data (default = `FALSE`)

- returnInfo:

  Return all the data used in SDA (default = `FALSE`)

- silent:

  Silent-ish mode (default = `FALSE`)

- noTitle:

  Do not generate a title (only the subtitle) (default = `FALSE`)

- tsName:

  Name of y added as a subtitle to the plot (default = `"y"`)

## Value

A list object containing:

- A data matrix `PLAW` with columns `freq.norm`, `size` and `bulk`.

- Estimate of scaling exponent `sap` based on a fit over the standard
  range (`fullRange`), or on a user defined range `fitRange`.

- Estimate of the the Fractal Dimension (`FD`) using conversion
  formula's reported in Hasselman(2013).

- Information output by various functions.

## References

Hasselman, F. (2013). When the blind curve is finite: dimension
estimation and model inference based on empirical waveforms. Frontiers
in Physiology, 4, 75. https://doi.org/10.3389/fphys.2013.00075

## See also

Other Fluctuation Analyses: [`fd_RR()`](fd_RR.md),
[`fd_allan()`](fd_allan.md), [`fd_dfa()`](fd_dfa.md),
[`fd_mfdfa()`](fd_mfdfa.md), [`fd_psd()`](fd_psd.md),
[`fd_sev()`](fd_sev.md)

## Author

Fred Hasselman
