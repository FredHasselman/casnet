# Symbolic representation

Return a discrete representation of `y` by transforming it into an
unordered categorical variable which indicates whether a value went up,
down, or remained the same relative to the previous value.

## Usage

``` r
ts_symbolic(y, keepNA = TRUE, usePlateaus = FALSE, doPlot = FALSE)
```

## Arguments

- y:

  Numeric vector or matrix to be discretised. Will return attributes
  with labels for each column in the matrix.

- keepNA:

  If `TRUE`, any `NA` values will first be removed and later re-inserted
  into the discretised time series. (default = `TRUE`)

- usePlateaus:

  Give consecutive `"same"` values after `"peak"` or `"trough"` a
  `"peak"` or `"trough"` label instrad of `"same"`. (default = `FALSE`)

- doPlot:

  Create a plot of the symbolized series. (default = `FALSE)`)

## Value

A symbolic version of `y`

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
[`ts_standardise()`](ts_standardise.md),
[`ts_sumorder()`](ts_sumorder.md), [`ts_trimfill()`](ts_trimfill.md),
[`ts_windower()`](ts_windower.md)
