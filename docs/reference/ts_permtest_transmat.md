# Permutation Test: Transition Matrix

Monte Carlo resampling of a time series using a discretised version of
`y`, a sequence of `bin` numbers with unique values equal to `nbins`:

1.  The discrete version of `y` will be used to generate a transition
    matrix of size `nbins X nbins`.

2.  This transition matrix will be used to resample values

## Usage

``` r
ts_permtest_transmat(
  y1,
  y2 = NULL,
  targetValue = 0,
  nbins = ceiling(2 * length(y1)^(1/3)),
  Nperms = 19,
  alpha = 0.05,
  keepNA = TRUE
)
```

## Arguments

- y1:

  Time series 1. The goal of the permutation test will be to decide
  whether the difference `y1-targetValue != 0` for each time point,
  given `alpha`.

- y2:

  An optional second time series. If this timeseries is provided then
  the goal of the permutation test will be the to decide wether the
  difference `y2-y1 != targetValue` for each time point, given `alpha`.

- targetValue:

  The target value for the permutation test. If `NULL`, the function
  will return a data frame with the block randomised surrogates columns
  (default = `0`)

- nbins:

  Number of bins to use (default = `ceiling(2*length(y1)^(1/3))`)

- Nperms:

  Number of permutations (default = `19`)

- alpha:

  Alpha level for deciding significance (default = `0.05`)

- keepNA:

  keepNA

## Value

Resampled series

## See also

Other Time series operations: [`ts_center()`](ts_center.md),
[`ts_changeindex()`](ts_changeindex.md),
[`ts_checkfix()`](ts_checkfix.md), [`ts_detrend()`](ts_detrend.md),
[`ts_diff()`](ts_diff.md), [`ts_discrete()`](ts_discrete.md),
[`ts_duration()`](ts_duration.md), [`ts_embed()`](ts_embed.md),
[`ts_integrate()`](ts_integrate.md), [`ts_levels()`](ts_levels.md),
[`ts_peaks()`](ts_peaks.md),
[`ts_permtest_block()`](ts_permtest_block.md),
[`ts_rasterize()`](ts_rasterize.md), [`ts_sd()`](ts_sd.md),
[`ts_slice()`](ts_slice.md), [`ts_slopes()`](ts_slopes.md),
[`ts_standardise()`](ts_standardise.md),
[`ts_sumorder()`](ts_sumorder.md), [`ts_symbolic()`](ts_symbolic.md),
[`ts_trimfill()`](ts_trimfill.md), [`ts_windower()`](ts_windower.md)

## Examples

``` r
set.seed(4321)
y <- rnorm(5000)
ts_permtest_transmat(y)
#>       rbin          ry
#>  [1,]   16  0.02075621
#>  [2,]   20  0.78987315
#>  [3,]   11 -1.14036837
#>  [4,]    9 -1.46663841
#>  [5,]   14 -0.38170474
#>  [6,]   17  0.13077913
#>  [7,]   22  1.21152759
#>  [8,]   11 -1.01566755
#>  [9,]   16  0.03722538
#> [10,]   17  0.28382814
#> [11,]   20  0.73169567
#> [12,]   12 -0.79037980
#> [13,]   13 -0.65595026
#> [14,]   12 -0.82335198
#> [15,]   19  0.67282503
#> [16,]   25  1.84039798
#> [17,]   15 -0.10940063
#> [18,]    7 -1.78496950
#> [19,]   19  0.70681715
#> [20,]   21  0.96253150
#> [21,]   22  1.31323045
#> [22,]   13 -0.69397775
#> [23,]   10 -1.22548324
#> [24,]   15 -0.25772863
#> [25,]   20  0.71951206
#> [26,]   21  0.93210270
#> [27,]   13 -0.66000699
#> [28,]    9 -1.41892332
#> [29,]   19  0.51425756
#> [30,]   24  1.73325662
#> [31,]   11 -0.99692179
#> [32,]   23  1.36824680
#> [33,]   18  0.33640892
#> [34,]   26  1.96289025
#> [35,]   26  2.01432372
```
