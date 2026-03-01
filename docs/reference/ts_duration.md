# Time series to Duration series

Time series to Duration series

## Usage

``` r
ts_duration(
  y,
  timeVec = stats::time(y),
  fs = stats::frequency(y),
  tolerance = 0
)
```

## Arguments

- y:

  A time series, numeric vector, or categorical variable.

- timeVec:

  A vector, same length as `y` containing timestamps, or, sample
  indices.

- fs:

  Optional sampling frequency if timeVec represents sample indices. An
  extra column `duration.fs` will be added which represents
  `1/fs * duration in samples`

- tolerance:

  A number `tol` indicating a range `[y-tol,y+tol]` to consider the same
  value. Useful when `y` is continuous (`default = 0`)

## Value

A data frame

## See also

Other Time series operations: [`ts_center()`](ts_center.md),
[`ts_changeindex()`](ts_changeindex.md),
[`ts_checkfix()`](ts_checkfix.md), [`ts_detrend()`](ts_detrend.md),
[`ts_diff()`](ts_diff.md), [`ts_discrete()`](ts_discrete.md),
[`ts_embed()`](ts_embed.md), [`ts_integrate()`](ts_integrate.md),
[`ts_levels()`](ts_levels.md), [`ts_peaks()`](ts_peaks.md),
[`ts_permtest_block()`](ts_permtest_block.md),
[`ts_permtest_transmat()`](ts_permtest_transmat.md),
[`ts_rasterize()`](ts_rasterize.md), [`ts_sd()`](ts_sd.md),
[`ts_slice()`](ts_slice.md), [`ts_slopes()`](ts_slopes.md),
[`ts_standardise()`](ts_standardise.md),
[`ts_sumorder()`](ts_sumorder.md), [`ts_symbolic()`](ts_symbolic.md),
[`ts_trimfill()`](ts_trimfill.md), [`ts_windower()`](ts_windower.md)

## Examples

``` r
library(invctr)
# Create data with events and their timecodes
coder <- data.frame(beh=c("stare","stare","coffee","type","type","stare"),t=c(0,5,10,15,20,25))

ts_duration(y = coder$beh, timeVec = coder$t)
#>   y y.name t.start t.end duration.time duration.samples duration.fs
#> 1 1  stare       0    10            10                2           2
#> 2 2 coffee      10    15             5                1           1
#> 3 3   type      15    25            10                2           2
#> 4 1  stare      25    25             0                1           1
```
