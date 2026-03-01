# Detect levels in time series

Use recursive partitioning function
[`rpart::rpart()`](https://rdrr.io/pkg/rpart/man/rpart.html) to perform
a 'classification' of relatively stable levels in a timeseries.

## Usage

``` r
ts_levels(
  y,
  minDataSplit = NROW(y)/5,
  minLevelDuration = round(minDataSplit/3),
  changeSensitivity = 0.01,
  maxLevels = 30,
  method = c("anova", "poisson", "class", "exp")[1],
  crossValidations = 10,
  minChange = NA,
  returnTrends = FALSE,
  Trend_minDataSplit = minDataSplit,
  Trend_minLevelDuration = minLevelDuration,
  Trend_changeSensitivity = changeSensitivity,
  Trend_maxLevels = maxLevels,
  Trend_method = method,
  Trend_crossValidations = crossValidations,
  Trend_minChange = NA,
  doLevelPlot = FALSE,
  doTreePlot = FALSE,
  returnTree = FALSE,
  returnPlot = FALSE,
  silent = FALSE
)
```

## Arguments

- y:

  A time series of numeric vector

- minDataSplit:

  An integer indicating how many datapoints should be in a segment
  before it will be analysed for presence of a level or trend change
  (default = `12`)

- minLevelDuration:

  Minimum duration (number of datapoint) of a level (default =
  `round(minDataSplit/3)`)

- changeSensitivity:

  A number indicating a criterion of change that must occur before
  declaring a new level. Higher numbers indicate higher levels of change
  must occur before a new level is considered. For example, if
  `method = "anova"`, the overall `R^2` after a level is introduced must
  increase by the value of `changeSensitivity`, see the `cp` parameter
  in
  [`rpart::rpart.control()`](https://rdrr.io/pkg/rpart/man/rpart.control.html).
  (default = `0.01`)

- maxLevels:

  Approximately the maximum number of levels tht will be detectd. The
  value indicates the node-depth of the final tree with the root node
  being depth 0, see
  [`rpart::rpart.control()`](https://rdrr.io/pkg/rpart/man/rpart.control.html)
  argument `maxdepth` (default = 30)

- method:

  The partitioning method to use, see the manual pages of
  [rpart::rpart](https://rdrr.io/pkg/rpart/man/rpart.html) for details.

- crossValidations:

  The number of cross-validations (default = `10`)

- minChange:

  After the call to
  [rpart::rpart](https://rdrr.io/pkg/rpart/man/rpart.html), adjust
  detected level changes to a minimum absolute change in `y`, e.g.
  `sd(y)` or in case of a discrete scale, a minimal scale change that
  represents significant level change relative to the interpretation of
  the scale. If a level change is smaller than `minChange`, the previous
  level will be continued. Note that this is an iterative process
  starting at the beginning of the series and 'correcting' towards the
  end. The results are stored in `p_adj`. Set to `NA` to skip, which
  means `p_adj` will be identical to `p` (default = `NA`)

- returnTrends:

  Should stationary trends also be estimated and returned? Unless
  otherwise specified, the arguments will be the same as for the level
  detection and argument `minChange` will be ignored for trends. This
  command will just run `ts_levels(diff(y))` and add the slope segments
  and values to the level plot. See `examples` on how to create a custom
  plot (default = `FALSE`)

- Trend_minDataSplit:

  see `minDataSplit`

- Trend_minLevelDuration:

  see `minLevelDuration`

- Trend_changeSensitivity:

  see `changeSensitivity`

- Trend_maxLevels:

  see `maxLevels`

- Trend_method:

  see `method`

- Trend_crossValidations:

  see `crossValidations`

- Trend_minChange:

  see `minChange`

- doLevelPlot:

  Should a plot with the original series and the levels and/or trends be
  produced? If `returnTrends = TRUE` sloped regions will be displayed as
  segments (default = `FALSE`)

- doTreePlot:

  Should a plot of the decision tree be produced. If
  `returnTrends = TRUE` there will be 2 trees. This requires package
  [partykit](https://cran.r-project.org/web/packages/partykit/index.html).
  Use [grid::grid.grab](https://rdrr.io/r/grid/grid.grab.html) to grab
  the tree plot object as a grob (default = `FALSE`)

- returnTree:

  should the tree object from
  [rpart::rpart](https://rdrr.io/pkg/rpart/man/rpart.html) be returned
  in the output? (default = `FALSE`)

- returnPlot:

  if `TRUE` returns the levelplot as a
  [ggplot2::ggplot2](https://ggplot2.tidyverse.org/reference/ggplot2-package.html)
  object (default = `FALSE`)

- silent:

  silent(-ish) mode

## Value

A list object with fields `tree` and `pred`. The latter is a data frame
with columns `x` (time), `y` (the variable of interest) and `p` the
predicted levels in `y` and `p_adj`, the levels in `p` but adjusted for
the value passed to `minChange`.

## See also

Other Time series operations: [`ts_center()`](ts_center.md),
[`ts_changeindex()`](ts_changeindex.md),
[`ts_checkfix()`](ts_checkfix.md), [`ts_detrend()`](ts_detrend.md),
[`ts_diff()`](ts_diff.md), [`ts_discrete()`](ts_discrete.md),
[`ts_duration()`](ts_duration.md), [`ts_embed()`](ts_embed.md),
[`ts_integrate()`](ts_integrate.md), [`ts_peaks()`](ts_peaks.md),
[`ts_permtest_block()`](ts_permtest_block.md),
[`ts_permtest_transmat()`](ts_permtest_transmat.md),
[`ts_rasterize()`](ts_rasterize.md), [`ts_sd()`](ts_sd.md),
[`ts_slice()`](ts_slice.md), [`ts_slopes()`](ts_slopes.md),
[`ts_standardise()`](ts_standardise.md),
[`ts_sumorder()`](ts_sumorder.md), [`ts_symbolic()`](ts_symbolic.md),
[`ts_trimfill()`](ts_trimfill.md), [`ts_windower()`](ts_windower.md)

## Author

Fred Hasselman

## Examples

``` r
# Levels in white noise?

set.seed(4321)
y <- rnorm(100)
wn <- ts_levels(y)
#> Skipping adjustment by argument minChange...
plot(wn$pred$x,wn$pred$y, type = "l")
lines(wn$pred$p, col = "red3", lwd = 2)

# This is due to the default changeSensitivity of 0.01

wn2 <- ts_levels(y,changeSensitivity = .1)
#> Skipping adjustment by argument minChange...
lines(wn2$pred$p, col = "steelblue", lwd = 2)


```
