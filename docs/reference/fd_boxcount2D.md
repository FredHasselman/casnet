# 2D Boxcount for 1D signal

2D Boxcount for 1D signal

## Usage

``` r
fd_boxcount2D(
  y = NA,
  unitSquare = TRUE,
  image2D = NA,
  resolution = 1,
  removeTrend = FALSE,
  polyOrder = 1,
  standardise = c("none", "mean.sd", "median.mad")[1],
  adjustSumOrder = FALSE,
  scaleMin = 0,
  scaleMax = floor(log2(NROW(y) * resolution)),
  scaleS = NA,
  dataMin = 2^(scaleMin + 1),
  maxData = 2^(scaleMax - 1),
  doPlot = FALSE,
  returnPlot = FALSE,
  returnPLAW = FALSE,
  returnInfo = FALSE,
  returnLocalScaling = FALSE,
  silent = FALSE,
  noTitle = FALSE,
  tsName = "y"
)
```

## Arguments

- y:

  A numeric vector or time series object.

- unitSquare:

  Create unit square image of `y`? This is required for estimating FD of
  time series (default = `TRUE`)

- image2D:

  A matrix representing a 2D image, argument `y` and `unitSquare` will
  be ignored (default = `NA`)

- resolution:

  The resolution used to embed the timeseries in 2D, a factor by which
  the dimensions the matrix will be multiplied (default = `1`)

- removeTrend:

  If `TRUE`, will call [ts_detrend](ts_detrend.md) on `y` (default =
  `FALSE`)

- polyOrder:

  Order of polynomial trend to remove if `removeTrend = `TRUE“

- standardise:

  Standardise `y` using [`ts_standardise()`](ts_standardise.md) with
  `adjustN = FALSE` (default = `none`)

- adjustSumOrder:

  Adjust the order of the time series (by summation or differencing),
  based on the global scaling exponent, see e.g.
  <https://www.frontiersin.org/files/Articles/23948/fphys-03-00141-r2/image_m/fphys-03-00141-t001.jpg>Ihlen (2012)
  (default = \`FALSE“)

- scaleMin:

  Minimium scale value (as `2^scale`) to use (default = `0`)

- scaleMax:

  Maximum scale value (as `2^scale`) to use (default = `max` of
  `log2(nrows)` and `log2(ncols)`)

- scaleS:

  If not `NA`, pass a numeric vector listing the scales (as a power of
  `2`) on which to evaluate the boxcount. Arguments `scaleMax`,
  `scaleMin`, and `scaleResolution` will be ignored (default = `NA`)

- maxData:

  Maximum number of time/data points inside a box for it to be included
  in the slope estimation (default = `2^scaleMax`)

- doPlot:

  Return the log-log scale versus bulk plot with linear fit (default =
  `TRUE`).

- returnPlot:

  Return ggplot2 object (default = `FALSE`)

- returnPLAW:

  Return the power law data (default = `FALSE`)

- returnInfo:

  Return all the data used in DFA (default = `FALSE`)

- returnLocalScaling:

  Return estimates of FD for each scale

- silent:

  Silent-ish mode (default = `TRUE`)

- noTitle:

  Do not generate a title (only the subtitle)

- tsName:

  Name of y added as a subtitle to the plot (default = `y`)

## Value

The boxcount fractal dimension and the 'local' boxcount fractal
dimension

## Note

This function was inspired by the `Matlab` function `boxcount.m`
[written by F.
Moisy](http://www.fast.u-psud.fr/~moisy/ml/boxcount/html/demo.md).
`Fred Hasselman` adapted the function for `R` for the purpose of the
unit square boxcount analysis for 1D time series. The original Matlab
toolbox has more options and contains more functions (e.g. `1D` and `3D`
boxcount).

## Examples

``` r
fd_boxcount2D(y = rnorm(100))
#> 
#> 
#> Raterizing time series... Done!
#> Performing 2D boxcount...Done!
#> 
#> ~~~o~~o~~casnet~~o~~o~~~
#> 
#>  2D boxcount of 1D curve 
#> 
#>  Full range (n = 9)
#> FD = 1.52 
#> 
#>  Exclude scales (n = 7)
#> FD = 1.52
#> 
#> ~~~o~~o~~casnet~~o~~o~~~

```
