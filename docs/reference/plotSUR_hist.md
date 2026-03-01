# Surrogate Test

Surrogate Test

## Usage

``` r
plotSUR_hist(
  surrogateValues,
  observedValue,
  sides = c("two.sided", "greater", "less")[1],
  binWidth = NULL,
  measureName = "",
  title = "",
  doPlot = TRUE,
  returnOnlyPvalue = FALSE
)
```

## Arguments

- surrogateValues:

  Vector of measures based on surrogate time series

- observedValue:

  The measure obtained from the observed value

- sides:

  Is this a 1 or 2-sided test (default = `1`)

- binWidth:

  The size of the histogram bins. The default is to look for the max.
  number of digits and set the width to `1/10^(Ndigits-1)`. If integers
  are detectec width will be set to 1.

- measureName:

  Label for x-axis

- title:

  A title for the plot

- doPlot:

  Plot a histogram of the distribution (default = `TRUE`)

- returnOnlyPvalue:

  Do not return the graph, just the point p-value (default = `FALSE`)

  alpha Significance threshold for the test. This value is currently
  calculated from the data as \\\frac{1}{rank}\*Nsides\\, setting it
  will not have an effect.

## Value

A point p-value for the observed value, and/or a histogram of the
distribution (`ggplot2` object).
