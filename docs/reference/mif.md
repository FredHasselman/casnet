# Mutual Information Function

Calculate the lagged mutual information fucntion within (auto-mif) or
between (cross-mif) time series, or, conditional on another time series
(conditional-cross-mif). Alternatively, calculate the total information
of a multivariate dataset for different lags.

## Usage

``` r
mif(
  y,
  lags = -10:10,
  nbins = ceiling(2 * NROW(y)^(1/3)),
  doPlot = FALSE,
  surTest = FALSE,
  alpha = 0.05
)
```

## Arguments

- y:

  A `Nx1` matrix for auto-mif, a `Nx2` matrix or data frame for
  cross-mif, a `Nx3` matrix or data frame for mif between col 1 and 2
  conditional on col 3; or a `NxM` matrix or data frame for the
  multi-information function. Mutual information for each lag will be
  calculated using functions in package
  [`infotheo::infotheo()`](https://rdrr.io/pkg/infotheo/man/infotheo.html)
  for `lags` lagged versions of the time series.

- lags:

  The lags to evaluate mutual information.

- nbins:

  The number of bins passed to
  [`infotheo::discretize()`](https://rdrr.io/pkg/infotheo/man/discretize.html)
  if y is a matrix or [`ts_discrete()`](ts_discrete.md)

- doPlot:

  Produce a plot of the symbolic time series by calling
  [`plotRED_mif()`](plotRED_mif.md) (default = `FALSE`)

- surTest:

  If `TRUE`, a surrogate will be conducted using simple surrogates. The
  surrogates will be created from the transition probabilities of the
  discretised time series, i.e. the probability of observing bin `j`
  when the current value is in bin `j`. The number of surrogates needed
  will be computed based on the value of the `alpha` parameter,
  conceived as a one-sided test: `mi > 0`.

- alpha:

  The alpha level for the surrogate test (default = `0.05`)

## Value

The auto- or cross-mi function

## See also

Other Redundancy measures (mutual information):
[`mi_interlayer()`](mi_interlayer.md), [`mi_mat()`](mi_mat.md)

## Examples

``` r
# Lags to evaluate mututal information
lags <- -10:30

# Auto-mutual information
y1 <- sin(seq(0, 100, by = 1/8)*pi)

(mif(data.frame(y1),lags = lags))
#>      -10       -9       -8       -7       -6       -5       -4       -3 
#> 1.656295 1.929039 1.645481 1.571829 1.567038 1.680948 1.460478 1.497499 
#>       -2       -1        0        1        2        3        4        5 
#> 1.476553 2.868931 2.868931 2.868931 1.476553 1.497499 1.460478 1.680948 
#>        6        7        8        9       10       11       12       13 
#> 1.567038 1.571829 1.645481 1.929039 1.656295 1.550729 1.924130 1.563809 
#>       14       15       16       17       18       19       20       21 
#> 1.607910 1.621651 1.560486 1.859941 1.474408 1.626296 1.476240 1.515003 
#>       22       23       24       25       26       27       28       29 
#> 1.453358 1.745546 1.566279 1.964918 1.604698 1.666703 1.567914 1.564823 
#>       30 
#> 1.734763 
#> attr(,"miType")
#> [1] "I(X;X)"
#> attr(,"lags")
#>  [1] -10  -9  -8  -7  -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6   7   8
#> [20]   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27
#> [39]  28  29  30
#> attr(,"nbins")
#> [1] 19

# Cross-mututal information, y2 is a lagged version y1
y2 <- y1[10:801]

y <- data.frame(ts_trimfill(y1, y2, action = "trim.cut"))
(mif(y,lags = lags))
#>      -10       -9       -8       -7       -6       -5       -4       -3 
#> 1.626296 1.474408 1.859941 1.560486 1.621651 1.607910 1.563809 1.924130 
#>       -2       -1        0        1        2        3        4        5 
#> 1.550729 1.656295 1.656295 1.656295 1.921974 1.687149 1.582642 1.610835 
#>        6        7        8        9       10       11       12       13 
#> 1.694766 1.457975 1.479514 1.467038 2.868171 1.483081 1.499478 1.474411 
#>       14       15       16       17       18       19       20       21 
#> 1.605462 1.578502 1.572016 1.660324 1.906882 1.644848 1.561052 1.923715 
#>       22       23       24       25       26       27       28       29 
#> 1.566571 1.604723 1.613996 1.696857 1.861352 1.475690 1.653503 1.467068 
#>       30 
#> 1.511635 
#> attr(,"miType")
#> [1] "I(X;Y)"
#> attr(,"lags")
#>  [1] -10  -9  -8  -7  -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6   7   8
#> [20]   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27
#> [39]  28  29  30
#> attr(,"nbins")
#> [1] 19

# Conditional mutual information, add some noise to y2 and add it as a 3rd column
y$s <- y2+rnorm(NROW(y2))
(mif(y,lags = lags))
#>      -10       -9       -8       -7       -6       -5       -4       -3 
#> 1.834929 1.736090 1.970503 1.833393 1.955936 1.964817 1.937624 2.039731 
#>       -2       -1        0        1        2        3        4        5 
#> 1.826218 1.881377 1.881377 1.881377 2.003972 1.933971 1.925013 1.982503 
#>        6        7        8        9       10       11       12       13 
#> 1.972004 1.826775 1.756128 1.733573 2.485812 1.819047 1.866779 1.838623 
#>       14       15       16       17       18       19       20       21 
#> 1.886196 1.910751 1.887256 1.857787 2.008641 1.966172 1.943974 2.106459 
#>       22       23       24       25       26       27       28       29 
#> 1.953648 1.914736 1.856459 1.881538 1.949608 1.809318 1.935573 1.857520 
#>       30 
#> 1.871952 
#> attr(,"miType")
#> [1] "I(X;Y|Z)"
#> attr(,"lags")
#>  [1] -10  -9  -8  -7  -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6   7   8
#> [20]   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27
#> [39]  28  29  30
#> attr(,"nbins")
#> [1] 19

# Multi-information, the information of the entire multivariate series at each lag
y$y3 <- cumsum(rnorm(NROW(y)))
(mif(y,lags = lags))
#>      -10       -9       -8       -7       -6       -5       -4       -3 
#> 4.993314 4.988131 4.985059 4.985194 5.130505 5.133889 5.125356 5.120978 
#>       -2       -1        0        1        2        3        4        5 
#> 5.120084 5.119594 5.119594 5.119594 5.120177 5.121040 5.125480 5.130071 
#>        6        7        8        9       10       11       12       13 
#> 5.131867 4.980106 4.981731 4.988115 4.991903 4.998527 4.996337 5.004096 
#>       14       15       16       17       18       19       20       21 
#> 5.007288 5.126603 5.125768 5.132465 5.135320 5.139722 5.140392 5.142012 
#>       22       23       24       25       26       27       28       29 
#> 5.146356 5.147014 5.146989 5.151418 4.997510 4.998788 5.009308 5.015651 
#>       30 
#> 5.011882 
#> attr(,"miType")
#> [1] "I(X;Y;Z;...;N)"
#> attr(,"lags")
#>  [1] -10  -9  -8  -7  -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6   7   8
#> [20]   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27
#> [39]  28  29  30
#> attr(,"nbins")
#> [1] 19

```
