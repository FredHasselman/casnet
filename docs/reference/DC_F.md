# Fluctuation Intensity

Fluctuation intensity is one of two components of which the product is
the Dynamic Complexity measure.

## Usage

``` r
dc_f(
  df,
  win = NROW(df),
  scale_min,
  scale_max,
  doPlot = FALSE,
  useVarNames = TRUE,
  colOrder = TRUE,
  useTimeVector = NA,
  timeStamp = "31-01-1999"
)
```

## Arguments

- df:

  A data frame containing multivariate time series data from 1 person.
  Rows should indicate time, columns should indicate the time series
  variables. All time series in `df` should be on the same scale, an
  error will be thrown if the range of the time series in`df` is not
  `[scale_min,scale_max]`.

- win:

  Size of window in which to calculate Dynamic Complexity. If
  `win < NROW(df)` the window will move along the time series with a
  stepsize of `1` (default = `NROW(df)`)

- scale_min:

  The theoretical minimum value of the scale. Used to calculate expected
  values, so it is important to set this to the correct value.

- scale_max:

  The theoretical maximum value of the scale. Used to calculate expected
  values, so it is important to set this to the correct value.

- doPlot:

  If `TRUE` shows a Complexity Resonance Diagram of the Dynamic
  Complexity and returns an invisible
  [`ggplot2::ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html)
  object. (default = `FALSE`)

- useVarNames:

  Use the column names of `df` as variable names in the Complexity
  Resonance Diagram (default = `TRUE`)

- colOrder:

  If `TRUE`, the order of the columns in `df` determines the of
  variables on the y-axis. Use `FALSE` for alphabetic/numeric order. Use
  `NA` to sort by by mean value of Dynamic Complexity (default = `TRUE`)

- useTimeVector:

  Parameter used for plotting. A vector of length `NROW(df)`, containing
  date/time information (default = `NA`)

- timeStamp:

  If `useTimeVector` is not `NA`, a character string that can be passed
  to
  [`lubridate::stamp()`](https://lubridate.tidyverse.org/reference/stamp.html)
  to format the the dates/times passed in `useTimeVector` (default =
  `"01-01-1999"`)

## Value

dataframe

## References

Haken H, & Schiepek G. (2006). *Synergetik in der Psychologie.
Selbstorganisation verstehen und gestalten*. Hogrefe, Göttingen.

Schiepek, G. (2003). A Dynamic Systems Approach to Clinical Case
Formulation. *European Journal of Psychological Assessment, 19*,
175-184. https://doi.org/10.1027//1015-5759.19.3.175

Schiepek, G., & Strunk, G. (2010). The identification of critical
fluctuations and phase transitions in short term and coarse-grained time
series-a method for the real-time monitoring of human change processes.
*Biological cybernetics, 102(3)*, 197-207.
https://doi.org/10.1007/s00422-009-0362-1

## See also

Use [`dc_win()`](dc_win.md) to get the dynamic complexity measure.

Other Dynamic Complexity functions: [`dc_ccp()`](dc_ccp.md),
[`dc_d()`](dc_d.md), [`dc_win()`](dc_win.md),
[`plotDC_ccp()`](plotDC_ccp.md), [`plotDC_lvl()`](plotDC_lvl.md),
[`plotDC_res()`](plotDC_res.md)

## Author

Merlijn Olthof

Fred Hasselman
