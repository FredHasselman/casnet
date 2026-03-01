# Detect slopes in time series

Use recursive partitioning function
[`rpart::rpart()`](https://rdrr.io/pkg/rpart/man/rpart.html) to perform
a 'classification' of relatively stable slopes in a timeseries.

## Usage

``` r
ts_slopes(
  y,
  minDataSplit = NROW(y)/4,
  minSlopeDuration = round(minDataSplit/3),
  changeSensitivity = 0.01,
  maxSlopes = floor(NROW(y)/minSlopeDuration),
  method = c("anova", "poisson", "class", "exp")[1],
  minChange = NA,
  doSlopePlot = FALSE,
  doTreePlot = FALSE
)
```

## Arguments

- y:

  A time series of numeric vector

- minDataSplit:

  An integer indicating how many datapoints should be in a segment
  before it will be analysed for presence of a slope (default = `12`)

- minSlopeDuration:

  Minimum duration (number of datapoints) of a slope (default =
  `round(minDataSplit/3)`)

- changeSensitivity:

  A number indicating a criterion of change that must occur before
  declaring the presence of a slope Higher numbers indicate higher
  levels of change must occur before a slope is considered. For example,
  if `method = "anova"`, the overall `R^2` after a slope is introduced
  must increase by the value of `changeSensitivity`, see the `cp`
  parameter in
  [`rpart::rpart.control()`](https://rdrr.io/pkg/rpart/man/rpart.control.html).
  (default = `0.01`)

- maxSlopes:

  Maximum number of levels in one series (default = floor(max(NROW(y),
  na.rm = TRUE)/minSlopeDuration))

- method:

  The partitioning method to use, see the manual pages of
  [rpart::rpart](https://rdrr.io/pkg/rpart/man/rpart.html) for details.

- minChange:

  After the call to
  [rpart::rpart](https://rdrr.io/pkg/rpart/man/rpart.html), adjust
  detected slope value to a minimum absolute change in `y`. If a slope
  value is smaller than `minChange`, the previous slope will be
  continued. Set e.g. to `sd(diff(y), na.rm = TRUE)`. Note that this is
  an iterative process starting at the beginning of the series and
  'correcting' towards the end. The results are stored in `p_adj`. Set
  to `NA` to skip, which means `p_adj` will be identical to `p` (default
  = `NA`)

- doSlopePlot:

  Should a plot with the original series and the levels be produced?
  (default = `FALSE`)

- doTreePlot:

  Should a plot of the decision tree be produced. This requires package
  [partykit](https://cran.r-project.org/web/packages/partykit/index.html)
  (default = `FALSE`)

## Value

A list object with fields `tree` and `pred`. The latter is a data frame
with columns `x` (time), `y` (the variable of interest) and `p` the
predicted slopes in `y` and `p_adj`, the slopes in `p` but adjusted for
the value passed to `minChange`.

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
[`ts_slice()`](ts_slice.md), [`ts_standardise()`](ts_standardise.md),
[`ts_sumorder()`](ts_sumorder.md), [`ts_symbolic()`](ts_symbolic.md),
[`ts_trimfill()`](ts_trimfill.md), [`ts_windower()`](ts_windower.md)

## Author

Fred Hasselman

## Examples

``` r
# Slopes in white noise?

set.seed(4321)
y <- rnorm(100)
wn <- ts_levels(y)
#> Skipping adjustment by argument minChange...
plot(wn$pred$x,wn$pred$y, type = "l")
lines(wn$pred$p, col = "red3", lwd = 2)

# This is due to the default changeSensitivity of 0.01

wn2 <- ts_slopes(y,changeSensitivity = .1)
#> Skipping adjustment by argument minChange...
lines(wn2$pred$p, col = "steelblue", lwd = 2)


```
