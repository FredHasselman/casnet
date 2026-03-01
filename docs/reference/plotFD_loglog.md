# Plot output from fluctuation analyses based on log-log regression

Plot output from fluctuation analyses based on log-log regression

## Usage

``` r
plotFD_loglog(
  fd.OUT,
  title = "",
  subtitle = "",
  xlabel = "Bin size",
  ylabel = "Fluctuation",
  logBase = NA,
  doPlot = TRUE,
  returnPlot = FALSE
)
```

## Arguments

- fd.OUT:

  Output from one of the `fd_` functions that use log-log regression to
  get scaling exponents.

- title:

  Plot title

- subtitle:

  Plot subtitle

- xlabel:

  x label

- ylabel:

  y label

- logBase:

  base of the log used

- doPlot:

  Display the plot (A plot object is always returned invisibly)

- returnPlot:

  return a
  [ggplot2::ggplot2](https://ggplot2.tidyverse.org/reference/ggplot2-package.html)
  object

## Value

A ggplot object
