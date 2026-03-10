# Sample Entropy

Sample Entropy

## Usage

``` r
inf_SampEn(
  y,
  m = 2,
  r = 0.2,
  D = NA,
  fs = NULL,
  standardise = c("none", "mean.sd", "median.mad")[1],
  transformBefore = TRUE,
  removeTrend = c("no", "poly", "adaptive", "bridge")[1],
  polyOrder = 1,
  relativeEntropy = FALSE,
  returnInfo = FALSE,
  silent = FALSE
)
```

## Arguments

- y:

  A numeric vector or time series object.

- m:

  The size of the window in which tho evaluate whether a pattern repeats
  (default = `2`)

- r:

  A factor that will determine the threshold for similarity of values,
  calculated as r x D (default = `0.2`)

- D:

  Commonly the standard deviation of the time series, the similarity
  threshold will be calculated as r x D. Note that if the series is
  detrended and/or standardised and `D = NA` the standard deviation will
  be calculated after the transformations (default = `NA`)

- fs:

  Sample rate

- standardise:

  Standardise the series using [`ts_standardise()`](ts_standardise.md)
  with `adjustN = FALSE` (default = "mean.sd")

- transformBefore:

  Detrend/standardise before coarse graining. If set to `FALSE`, each
  coarsegrained series will be detrended/standardised separately
  (default = `TRUE`)

- removeTrend:

  Method to use for global detrending (default = `"poly"`)

- polyOrder:

  Order of global polynomial trend to remove if `removeTrend = "poly"`.
  If `removeTrend = "adaptive"` polynomials `1` to `polyOrder` will be
  evaluated and the best fitting curve (R squared) will be removed
  (default = `1`)

- relativeEntropy:

  The relative entropy, SampEn / (-1 \* log(1/length(y))) will be
  returned (default = `FALSE`)

- returnInfo:

  Return all the data used in SDA (default = `FALSE`)

- silent:

  Silent-ish mode (default = `FALSE`)

## Value

The sample entropy (SampEn) of the time series y.

## See also

info_MSE

Other Information based complexity measures: [`inf_MSE()`](inf_MSE.md)

## Examples

``` r

y <- rnorm(100)

# Similarity threshold is r x D = 0.2 * sd(y)
inf_SampEn(y)
#> 
#> 
#> fd_SampEn:   Sample rate was set to 1.
#> 
#> [1] 2.433613
#> attr(,"Relative Entropy")
#> [1] 0.5284524

# Similarity threshold is r = 0.2
inf_SampEn(y, D = 1)
#> 
#> 
#> fd_SampEn:   Sample rate was set to 1.
#> 
#> [1] 2.463853
#> attr(,"Relative Entropy")
#> [1] 0.5350189
```
