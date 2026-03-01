# Trim or Fill Vectors

Trim the largest vector by cutting it, or filling it with `NA`. Fill the
shortest vector with padding.

## Usage

``` r
ts_trimfill(
  x,
  y,
  action = c("fill", "trim.cut", "trim.NA")[1],
  type = c("end", "center", "front")[1],
  padding = 0,
  silent = TRUE
)
```

## Arguments

- x:

  A numeric vector

- y:

  A numeric vector

- action:

  Use `"fill"` to fill the shortest vector with `padding` (default);
  `"trim.cut"` to trim the longest vector to the length of the shortest;
  `"trim.NA"` to fill the longest vector with `NA`. This is a shortcut
  for running `action = "trim.cut"` with `padding=NA`, which can be
  useful if one wants to match the shortest series, but preserve the
  original length of largest vector.

- type:

  Should trimming or filling take place at the `"end"` (default), or
  `"front"` of the vector? The option `"center"` will try to distribute
  trimming by `NA` or filling by `padding` evenly across the front and
  end of the vector.

- padding:

  A value to use for padding (default = `0`)

- silent:

  Run silent-ish

## Value

A list with two vectors of equal length.

## See also

il_mi

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
[`ts_sumorder()`](ts_sumorder.md), [`ts_symbolic()`](ts_symbolic.md),
[`ts_windower()`](ts_windower.md)

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
[`ts_sumorder()`](ts_sumorder.md), [`ts_symbolic()`](ts_symbolic.md),
[`ts_windower()`](ts_windower.md)

## Author

Fred Hasselman
