# Standardise a vector

Standardise a vector

## Usage

``` r
ts_standardise(
  y,
  na.rm = TRUE,
  keepNAvalues = TRUE,
  type = c("mean.sd", "median.mad")[1],
  adjustN = TRUE
)
```

## Arguments

- y:

  A time series or numeric vector

- na.rm:

  Set the `na.rm` field

- keepNAvalues:

  If `na.rm = TRUE` and `keepNAvalues = TRUE`, any `NA` values in `y`
  will be re-inserted after transformation.

- type:

  Center on the `"mean"` and divide by `sd` (default), or center on
  `"median"` and divide by `mad`

- adjustN:

  If `TRUE`, apply Bessel's correction (divide by `N-1`) or return the
  unadjusted `SD` (divide by `N`) (default = `TRUE`)

## Value

A standardised vector

## See also

Other Time series operations: [`ts_center()`](ts_center.md),
[`ts_changeindex()`](ts_changeindex.md),
[`ts_checkfix()`](ts_checkfix.md), [`ts_detrend()`](ts_detrend.md),
[`ts_diff()`](ts_diff.md), [`ts_discrete()`](ts_discrete.md),
[`ts_duration()`](ts_duration.md), [`ts_embed()`](ts_embed.md),
[`ts_integrate()`](ts_integrate.md), [`ts_levels()`](ts_levels.md),
[`ts_peaks()`](ts_peaks.md),
[`ts_permtest_block()`](ts_permtest_block.md),
[`ts_permtest_transmat()`](ts_permtest_transmat.md),
[`ts_rasterize()`](ts_rasterize.md), [`ts_sd()`](ts_sd.md),
[`ts_slice()`](ts_slice.md), [`ts_slopes()`](ts_slopes.md),
[`ts_sumorder()`](ts_sumorder.md), [`ts_symbolic()`](ts_symbolic.md),
[`ts_trimfill()`](ts_trimfill.md), [`ts_windower()`](ts_windower.md)

## Author

Fred Hasselman
