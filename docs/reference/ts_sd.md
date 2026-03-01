# Standard Deviation estimates

Calculates the population estimate of the standard deviation, or the
unadjusted standard deviation.

## Usage

``` r
ts_sd(y, na.rm = TRUE, type = c("Bessel", "unadjusted")[1], silent = TRUE)
```

## Arguments

- y:

  Time series or numeric vector

- na.rm:

  Remove missing values before calculations

- type:

  Apply Bessel's correction (divide by N-1) or return unadjusted SD
  (divide by N)

- silent:

  Silent-ish mode (default = `TRUE`)

## Value

Standard deviation of `y`

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
[`ts_rasterize()`](ts_rasterize.md), [`ts_slice()`](ts_slice.md),
[`ts_slopes()`](ts_slopes.md), [`ts_standardise()`](ts_standardise.md),
[`ts_sumorder()`](ts_sumorder.md), [`ts_symbolic()`](ts_symbolic.md),
[`ts_trimfill()`](ts_trimfill.md), [`ts_windower()`](ts_windower.md)
