# Detrend a time series

Detrend a time series

## Usage

``` r
ts_detrend(y, polyOrder = 1, adaptive = FALSE)
```

## Arguments

- y:

  A time series of numeric vector

- polyOrder:

  order Order of polynomial trend to remove

## Value

Residuals after detrending polynomial of order `order`

## See also

Other Time series operations: [`ts_center()`](ts_center.md),
[`ts_changeindex()`](ts_changeindex.md),
[`ts_checkfix()`](ts_checkfix.md), [`ts_diff()`](ts_diff.md),
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

## Author

Fred Hasselman
