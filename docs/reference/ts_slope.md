# Calculate Kendall's tau in sliding window or around change point.

Kendall's tau at different lags or change point ranges

## Usage

``` r
ts_slope(y, win = NA, changepoint = NA, doPlot = FALSE)
```

## Arguments

- y:

  A numeric vector

- win:

  Numeric vector with 1 or 2 values. If `changepoint != NA`, Kendall's
  tau will be calculated in a window around the change point. If one
  value is provided a symmetric window, otherwise `changepoint - win[1]`
  and `changepoint + win[2]`. If `changepoint == NA` Kendall's tau will
  be calculated in a sliding window.

- changepoint:

  If not `NA`, it has to be an index of `y`. If `win` is `NA` Kendall's
  tau will be calculated in `1:changepoint`, otherwise the values of
  `win` will be used to create a window. If both `changepoint` and `win`
  are `NA` Kendall's tau will be calculated on all of `y`.

- doPlot:

  provide a plot

## Value

A data frame

## Examples

``` r
ts_slope(y=rnorm(100), doPlot = TRUE)
#> Error in ts_slope(y = rnorm(100), doPlot = TRUE): object 'X' not found
```
