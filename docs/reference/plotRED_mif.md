# Plot various MI functions

Plot various MI functions

## Usage

``` r
plotRED_mif(
  mif.OUT = NULL,
  lags = attr(mif.OUT, "lags"),
  nbins = attr(mif.OUT, "nbins"),
  surTest = FALSE,
  alpha = 0.05,
  doPlot = TRUE,
  returnMIFun = TRUE
)
```

## Arguments

- mif.OUT:

  Output from function [`mif()`](mif.md)

- lags:

  The lags to evaluate mutual information.

- nbins:

  The number of bins passed to
  [`infotheo::discretize()`](https://rdrr.io/pkg/infotheo/man/discretize.html)
  if y is a matrix or [`ts_discrete()`](ts_discrete.md)

- surTest:

  If `TRUE`, a surrogate will be conducted using simple surrogates. The
  surrogates will be created from the transition probabilities of the
  discretised time series, i.e. the probability of observing bin `j`
  when the current value is in bin `j`. The number of surrogates needed
  will be computed based on the value of the `alpha` parameter,
  conceived as a one-sided test: `mi > 0`.

- alpha:

  The alpha level for the surrogate test (default = `0.05`)

- doPlot:

  Produce a plot of the symbolic time series by calling `plotRED_mif()`
  (default = `FALSE`)

- returnMIFun:

  Return the data

## Value

Either an invisible ggplot2 object r a list containing the plot and the
data

## See also

Other Plot redundancy functions: [`plotRED_acf()`](plotRED_acf.md)
