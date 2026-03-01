# Plot ACF and PACF

Plot ACF and PACF

## Usage

``` r
plotRED_acf(
  y,
  Lmax = max(round(NROW(y)/4), 10),
  alpha = 0.05,
  doPlot = TRUE,
  returnCorFun = FALSE
)
```

## Arguments

- y:

  A time series or numeric vector

- Lmax:

  Maximum number of lags

- alpha:

  Significance level

- doPlot:

  Plot output

- returnCorFun:

  Return the data

## Value

Either an invisible ggplot2 object r a list containing the plot and the
data

## See also

Other Plot redundancy functions: [`plotRED_mif()`](plotRED_mif.md)
