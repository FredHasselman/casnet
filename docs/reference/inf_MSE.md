# Multi-Scale Entropy

Calculate the Multi-Scale (Sample) Entropy of a time series.

## Usage

``` r
inf_MSE(
  y,
  summaryFunction = c("mean", "median", "min", "max")[1],
  retainLength = FALSE,
  m = 2,
  r = 0.2,
  D = NA,
  fs = NULL,
  transformBefore = TRUE,
  removeTrend = c("no", "poly", "adaptive", "bridge")[1],
  polyOrder = 1,
  standardise = c("none", "mean.sd", "median.mad")[1],
  adjustSumOrder = FALSE,
  scaleMin = 1,
  scaleMax = floor(NROW(y)/10),
  scaleS = NA,
  overlap = 0,
  relativeEntropy = FALSE,
  doPlot = FALSE,
  returnPlot = FALSE,
  returnPLAW = FALSE,
  returnInfo = FALSE,
  silent = FALSE,
  noTitle = FALSE,
  tsName = "y"
)
```

## Arguments

- y:

  A numeric vector

- summaryFunction:

  How should the data be summarized in the bins? Can be `"mean"`,
  `"median"`, `"min"`, `"max"`, or, `"maxFreq"`. Value `"maxFreq"` is
  for categorical data and will pick the most frequently occurring
  category within the bin (default = `"mean"`)

- retainLength:

  Return only the bin values (`FALSE`), or retain the length of the
  original series? (default = `TRUE`)

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

- standardise:

  Standardise the series using [`ts_standardise()`](ts_standardise.md)
  with `adjustN = FALSE` (default = "mean.sd")

- adjustSumOrder:

  Adjust the time series (summation or difference), based on the global
  scaling exponent, see e.g. [Ihlen
  (2012)](https://www.frontiersin.org/files/Articles/23948/fphys-03-00141-r2/image_m/fphys-03-00141-t001.jpg)
  (default = `FALSE`)

- scaleMin:

  Minimum scale (in data points) to use for log-log regression (default
  = `4`)

- scaleMax:

  Maximum scale (in data points) to use for log-log regression (default
  = `stats::nextn(floor(NROW(y)/4), factors = 2)`)

- scaleS:

  If not `NA`, it should be a numeric vector listing the scales on which
  to evaluate the detrended fluctuations. Arguments
  `scaleMax, scaleMin, scaleResolution` will be ignored (default = `NA`)

- overlap:

  A number in `[0 ... 1]` representing the amount of 'bin overlap' when
  calculating the fluctuation. This reduces impact of arbitrary time
  series begin and end points. If `length(y) = 1024` and overlap is
  `.5`, a scale of `4` will be considered a sliding window of size `4`
  with step-size `floor(.5 * 4) = 2`, so for scale `128` step-size will
  be `64` (default = `NA`)

- relativeEntropy:

  The relative entropy, SampEn / (-1 \* log(1/length(y))) will be
  returned (default = `FALSE`)

- doPlot:

  Output the log-log scale versus fluctuation plot with linear fit by
  calling function [`plotFD_loglog()`](plotFD_loglog.md) (default =
  `TRUE`)

- returnPlot:

  Return ggplot2 object (default = `FALSE`)

- returnPLAW:

  Return the power law data (default = `FALSE`)

- returnInfo:

  Return all the data used in SDA (default = `FALSE`)

- silent:

  Silent-ish mode (default = `FALSE`)

- noTitle:

  Do not generate a title (only the subtitle) (default = `FALSE`)

- tsName:

  Name of y added as a subtitle to the plot (default = `"y"`)

## Value

MSE

## Note

The transformation settings (detrending and/or normalisation),
`transformBefore` will determine the value of the product r x D if
`D = NA`.

## See also

Other Information based complexity measures:
[`inf_SampEn()`](inf_SampEn.md)

## Examples

``` r
y <- rnorm(100)

out <- inf_MSE(y, r=.5, scaleMin=1, scaleMax=10, returnInfo = TRUE, doPlot = FALSE)
#> 
#> ~~~o~~o~~casnet~~o~~o~~~
#> 
#>  Multi Scale Entropy 
#> 
#>  Full range (n = 10)
#> Slope = -0.14 
#> 
#>  Fit range (n = 7)
#> Slope = -0.12
#> 
#> ~~~o~~o~~casnet~~o~~o~~~

```
