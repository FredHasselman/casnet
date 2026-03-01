# Cumulative Complexity Peaks (CCP)

Computes significant peaks in the dynamic complexity time series.

## Usage

``` r
dc_ccp(
  df_win,
  alpha_item = 0.05,
  alpha_time = 0.05,
  doPlot = FALSE,
  useVarNames = TRUE,
  colOrder = TRUE,
  useTimeVector = NA,
  timeStamp = "31-01-1999",
  markID = NA,
  markIDcolour = "grey",
  markIDlabel = "Time points of interest marked grey",
  markIDalpha = 0.5,
  NAdates = 1:(win - 1),
  trimFirstWin = TRUE
)
```

## Arguments

- df_win:

  A data frame containing series of Dynamic Complexity values obtained
  by running function [`dc_win()`](dc_win.md)

- alpha_item:

  The significance level of the one-sided Z-test used to determine which
  peaks are `> 0`.

- alpha_time:

  The significance level of the one-sided Z-test used to determine if
  the number of significant peaks (as determined by `alpha_item`) at a
  specific time stamp are `> 0`.

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

- markID:

  Numeric vector of integers in the range
  `[length of window, length of timeseries]`. Vertical lines will be
  drawn at these indices (default = `NA`)

- markIDcolour:

  Colour of time point markers (default = `"red"`)

- markIDlabel:

  Label added to subtitle explaining time point markers (default =
  `Time points of interest marked red`)

- markIDalpha:

  Alpha of time point marker colour (default = `.5`)

- NAdates:

  Should some dates be considered `NA`? Provide a numerical vector with
  indices, the default is to set `1:(win-1)` to NA. (default =
  `1:(win-1)`)

- trimFirstWin:

  Display the first empty window (`1:win-1`)? (default = `TRUE`)

## Value

A list with a dataframe of binary complexity peak indices and a
cumulative complexity peak index, a CCP diagram.

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

Other Dynamic Complexity functions: [`dc_d()`](dc_d.md),
[`dc_f()`](dc_f.md), [`dc_win()`](dc_win.md),
[`plotDC_ccp()`](plotDC_ccp.md), [`plotDC_lvl()`](plotDC_lvl.md),
[`plotDC_res()`](plotDC_res.md)

## Author

Merlijn Olthof

Fred Hasselman
