# Delay embedding of a time series

Create a state vector based on an embedding lag and a number of
embedding dimanesions.

## Usage

``` r
ts_embed(
  y,
  emDim,
  emLag,
  returnOnlyIndices = FALSE,
  doEmbed = TRUE,
  silent = TRUE
)
```

## Arguments

- y:

  Time series

- emDim:

  Embedding dimension

- emLag:

  Embedding lag

- returnOnlyIndices:

  Return only the index of y for each surrogate dimension, not the
  values (default = `FALSE`)

- doEmbed:

  Should the series be embedded? If `FALSE` adds attributes.

- silent:

  Silent-ish mode

## Value

The lag embedded time series

## Note

If `emLag = 0`, the assumption is the columns in `y` represent the
dimensions and `y` will be returned with attributes `emLag = 0` and
`emDim = NCOL(y)`. If `emLag > 0` and `NCOL(y)>1` the first column of
`y` will used for embedding and a warning will be triggered.

## See also

Other Time series operations: [`ts_center()`](ts_center.md),
[`ts_changeindex()`](ts_changeindex.md),
[`ts_checkfix()`](ts_checkfix.md), [`ts_detrend()`](ts_detrend.md),
[`ts_diff()`](ts_diff.md), [`ts_discrete()`](ts_discrete.md),
[`ts_duration()`](ts_duration.md), [`ts_integrate()`](ts_integrate.md),
[`ts_levels()`](ts_levels.md), [`ts_peaks()`](ts_peaks.md),
[`ts_permtest_block()`](ts_permtest_block.md),
[`ts_permtest_transmat()`](ts_permtest_transmat.md),
[`ts_rasterize()`](ts_rasterize.md), [`ts_sd()`](ts_sd.md),
[`ts_slice()`](ts_slice.md), [`ts_slopes()`](ts_slopes.md),
[`ts_standardise()`](ts_standardise.md),
[`ts_sumorder()`](ts_sumorder.md), [`ts_symbolic()`](ts_symbolic.md),
[`ts_trimfill()`](ts_trimfill.md), [`ts_windower()`](ts_windower.md)

Other Time series operations: [`ts_center()`](ts_center.md),
[`ts_changeindex()`](ts_changeindex.md),
[`ts_checkfix()`](ts_checkfix.md), [`ts_detrend()`](ts_detrend.md),
[`ts_diff()`](ts_diff.md), [`ts_discrete()`](ts_discrete.md),
[`ts_duration()`](ts_duration.md), [`ts_integrate()`](ts_integrate.md),
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
