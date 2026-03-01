# Slice a Matrix

Slices rows of a matrix into a list of matrices representing epochs of
length `epochSz`.

## Usage

``` r
ts_slice(y, epochSz = 4, overlap = NA, removeUnequal = FALSE)
```

## Arguments

- y:

  A matrix with timeseries as columns

- epochSz:

  Epoch size

- removeUnequal:

  Do not return bins whose length is not equal to `epochSz` (default =
  `FALSE`)

## Value

A list with epochs

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
[`ts_slopes()`](ts_slopes.md), [`ts_standardise()`](ts_standardise.md),
[`ts_sumorder()`](ts_sumorder.md), [`ts_symbolic()`](ts_symbolic.md),
[`ts_trimfill()`](ts_trimfill.md), [`ts_windower()`](ts_windower.md)

## Author

Fred Hasselman
