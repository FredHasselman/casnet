# Check and/or Fix a vector

Check and/or Fix a vector

## Usage

``` r
ts_checkfix(
  y,
  checkNumericVector = TRUE,
  checkWholeNumbers = FALSE,
  checkTimeVector = FALSE,
  checkPow2 = FALSE,
  checkScale = FALSE,
  checkSummationOrder = FALSE,
  checkNonStationarity = FALSE,
  checkNonHomogeneity = FALSE,
  fixNumericVector = FALSE,
  fixWholeNumbers = FALSE,
  fixTimeVector = FALSE,
  fixPow2 = FALSE,
  fixNA = TRUE,
  fixScale = FALSE,
  fixSummationOrder = FALSE,
  fixNonStationarity = FALSE,
  fixNonHomogeneity = FALSE
)
```

## Arguments

- y:

  A time series object or numeric vector

- checkNumericVector:

  is 1D numeric vector?

- checkWholeNumbers:

  contains only wholenumbers?

- checkTimeVector:

  has time vector?

- checkPow2:

  length is power of 2?

- checkScale:

  checkScale

- checkSummationOrder:

  checkSummationOrder

- checkNonStationarity:

  checkNonStationarity

- checkNonHomogeneity:

  checkNonHomogeneity

- fixNumericVector:

  return a 1D numeric vector (WARNING: Data frames and Matrices with
  NCOL \> 1 wil be converted to long form)

- fixWholeNumbers:

  fixWholeNumber

- fixTimeVector:

  fixTimeVector

- fixPow2:

  foxPow2

- fixNA:

  fixNA

- fixScale:

  fixScale

- fixSummationOrder:

  fixSummationOrder

- fixNonStationarity:

  fixNonStationarity

- fixNonHomogeneity:

  fixNonHomogeneity

## Value

A 'check' report and/or a 'fixed' vector y.

## See also

Other Time series operations: [`ts_center()`](ts_center.md),
[`ts_changeindex()`](ts_changeindex.md),
[`ts_detrend()`](ts_detrend.md), [`ts_diff()`](ts_diff.md),
[`ts_discrete()`](ts_discrete.md), [`ts_duration()`](ts_duration.md),
[`ts_embed()`](ts_embed.md), [`ts_integrate()`](ts_integrate.md),
[`ts_levels()`](ts_levels.md), [`ts_peaks()`](ts_peaks.md),
[`ts_permtest_block()`](ts_permtest_block.md),
[`ts_permtest_transmat()`](ts_permtest_transmat.md),
[`ts_rasterize()`](ts_rasterize.md), [`ts_sd()`](ts_sd.md),
[`ts_slice()`](ts_slice.md), [`ts_slopes()`](ts_slopes.md),
[`ts_standardise()`](ts_standardise.md),
[`ts_sumorder()`](ts_sumorder.md), [`ts_symbolic()`](ts_symbolic.md),
[`ts_trimfill()`](ts_trimfill.md), [`ts_windower()`](ts_windower.md)
