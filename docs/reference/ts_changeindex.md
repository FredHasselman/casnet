# Find change indices

Find change indices

## Usage

``` r
ts_changeindex(
  y,
  returnRectdata = FALSE,
  groupVar = NULL,
  labelVar = NULL,
  discretize = FALSE,
  nbins = 5
)
```

## Arguments

- y:

  An indicator variable representing different levels of a variable or
  factor

- returnRectdata:

  Return a dataframe suitable for shading a `ggplot2` graph with
  [`ggplot2::geom_rect()`](https://ggplot2.tidyverse.org/reference/geom_tile.html)

- groupVar:

  Pass a value (length 1) or variable (length of y) that can be used as
  a variable to join the indices by if `returnRectdata = TRUE`

- labelVar:

  If `y` is not a character vector, provide a vector of labels equal to
  `length(y)`

- discretize:

  If `y` is a continuous variable, setting `discretize = TRUE` will
  partition the values of `y` into `nbins` number of bins, each value of
  `y` will be replaced by its bin number.

- nbins:

  Number of bins to use to change a continuous `y` (if
  `discretize = TRUE`) into a variable with `nbins` levels

## Value

Either a vector with the indices of change in `y`, or, a data frame with
variables `xmin,xmax,ymin,ymax,label`

## See also

Other Time series operations: [`ts_center()`](ts_center.md),
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
[`ts_sumorder()`](ts_sumorder.md), [`ts_symbolic()`](ts_symbolic.md),
[`ts_trimfill()`](ts_trimfill.md), [`ts_windower()`](ts_windower.md)

## Examples

``` r
 library(ggplot2)

 set.seed(1234)
 yy     <- noise_powerlaw(standardise = TRUE, N=50, alpha = -1)
 tr     <- ts_levels(yy, doTreePlot = TRUE)
#> Skipping adjustment by argument minChange...

 breaks <- ts_changeindex(tr$pred$p, returnRectdata = TRUE)
 breaks$cols <- casnet::getColours(length(breaks$label))

 ggplot(tr$pred) +
   geom_rect(data = breaks,
   aes(xmin = xmin, xmax=xmax, ymin=ymin, ymax=ymax, colour = label, fill = label),
   alpha = .3) +
   scale_colour_manual(values = breaks$cols) +
   scale_fill_manual(values = breaks$cols) +
   scale_x_continuous("time", expand = c(0,0)) +
   geom_line(aes(x=x,y=y)) +
   geom_step(aes(x=x,y=p), colour = "red3", size=1) +
   theme_bw() + theme(panel.grid.minor = element_blank())

```
