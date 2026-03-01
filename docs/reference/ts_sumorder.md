# Adjust time series by summation order

Many fluctuation analyses assume a time series' Hurst exponent is within
the range of `0.2 - 1.2`. If this is not the case it is sensible to make
adjustments to the time series, as well as the resutling Hurst exponent.

## Usage

``` r
ts_sumorder(y, scaleS = NULL, polyOrder = 1, dataMin = 4)
```

## Arguments

- y:

  A time series of numeric vector

- scaleS:

  The scales to consider for `DFA1`

- polyOrder:

  Order of polynomial for detrending in DFA (default = `1`)

- dataMin:

  Minimum number of data points in a bin needed to calculate detrended
  fluctuation

## Value

The input vector, possibly adjusted based on `H` with an attribute
`"Hadj"` containing an integer by which a Hurst exponent calculated from
the series should be adjusted.

## Details

Following recommendations by
<https://www.frontiersin.org/files/Articles/23948/fphys-03-00141-r2/image_m/fphys-03-00141-t001.jpg>Ihlen
(2012), a global Hurst exponent is estimated using DFA and `y` is
adjusted accordingly:

- `1.2 < H < 1.8` first derivative of y, atribute `Hadj = 1`

- `H > 1.8` second derivative of y, atribute `Hadj = 2`

- `H < 0.2` y is centered and integrated, atribute `Hadj = -1`

- `0.2 <= H <= 1.2 ` y is unaltered, atribute `Hadj = 0`

## References

Ihlen, E. A. F. E. (2012). Introduction to multifractal detrended
fluctuation analysis in Matlab. Frontiers in physiology, 3, 141.

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
[`ts_symbolic()`](ts_symbolic.md), [`ts_trimfill()`](ts_trimfill.md),
[`ts_windower()`](ts_windower.md)
