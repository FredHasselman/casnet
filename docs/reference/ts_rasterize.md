# Turn a 1D time series vector into a 2D curve

Turn a 1D time series vector into a 2D curve

## Usage

``` r
ts_rasterize(y, unitSquare = FALSE, toSparse = TRUE, resolution = 2)
```

## Arguments

- y:

  A 1D time series object or numeric vector.

- unitSquare:

  Convert the series to a unit square? (default = `FALSE`)

- toSparse:

  Convert to sparse Matrix (default = `FALSE`)

- resolution:

  Factor by which dimensions will be multiplied (default = `2`)

## Value

A (sparse) matrix representing the time series as a curve in 2D space

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
[`ts_sd()`](ts_sd.md), [`ts_slice()`](ts_slice.md),
[`ts_slopes()`](ts_slopes.md), [`ts_standardise()`](ts_standardise.md),
[`ts_sumorder()`](ts_sumorder.md), [`ts_symbolic()`](ts_symbolic.md),
[`ts_trimfill()`](ts_trimfill.md), [`ts_windower()`](ts_windower.md)

## Examples

``` r
# \donttest{
y <- rnorm(100)
plot(ts(y))


y_img <- ts_rasterize(y)
image(y_img,col=c("white","black"))# }


```
