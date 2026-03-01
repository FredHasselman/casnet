# Plot Complexity Resonance Diagram

Plot Complexity Resonance Diagram

## Usage

``` r
plotDC_res(
  df_win,
  win,
  useVarNames = TRUE,
  colOrder = TRUE,
  useTimeVector = NA,
  timeStamp = "31-01-1999",
  markID = NA,
  markIDcolour = "grey",
  markIDlabel = "Time points of interest marked grey",
  markIDalpha = 0.5,
  doPlot = TRUE,
  PlotMeanDC = TRUE,
  title = "Complexity Resonance Diagram",
  resVariable = "Dynamic Complexity",
  subtitle = "",
  xlabel = "Time",
  ylabel = "",
  NAdates = 1:win,
  trimFirstWin = TRUE
)
```

## Arguments

- df_win:

  A data frame containing series of Dynamic Complexity values obtained
  by running function [`dc_win()`](dc_win.md)

- win:

  Size of window in which to calculate Dynamic Complexity. If
  `win < NROW(df)` the window will move along the time series with a
  stepsize of `1` (default = `NROW(df)`)

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

- doPlot:

  If `TRUE` shows a Complexity Resonance Diagram of the Dynamic
  Complexity and returns an invisible
  [`ggplot2::ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html)
  object. (default = `FALSE`)

- title:

  A title for the plot.

- resVariable:

  Variable displayed in the plot.

- subtitle:

  A subtitle for the plot.

- xlabel:

  A label for the x-axis.

- ylabel:

  A label for the y-axis.

- NAdates:

  Should some dates be considered `NA`? Provide a numerical vector with
  indices, the default is to set `1:(win-1)` to NA. (default =
  `1:(win-1)`)

- trimFirstWin:

  Display the first empty window (`1:win-1`)? (default = `TRUE`)

## Value

An invisible ggplot2 object.

## See also

Other Dynamic Complexity functions: [`dc_ccp()`](dc_ccp.md),
[`dc_d()`](dc_d.md), [`dc_f()`](dc_f.md), [`dc_win()`](dc_win.md),
[`plotDC_ccp()`](plotDC_ccp.md), [`plotDC_lvl()`](plotDC_lvl.md)
