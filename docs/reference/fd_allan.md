# Allan Variance Analysis

Allan Variance Analysis

## Usage

``` r
fd_allan(
  y,
  fs = stats::tsp(stats::hasTsp(y))[3],
  useSD = FALSE,
  doPlot = FALSE,
  returnPlot = FALSE,
  returnPLAW = FALSE,
  returnInfo = FALSE,
  silent = FALSE,
  noTitle = FALSE,
  tsName = "y"
)
```

## Arguments

- y:

  A numeric vector or time series object

- fs:

  Sample frequency in Hz

- useSD:

  Use the standard deviation instead of variance?

- doPlot:

  Return the log-log scale versus fluctuation plot with linear fit
  (default = `TRUE`).

- returnPlot:

  Return ggplot2 object (default = `FALSE`)

- returnPLAW:

  Return the power law data (default = `FALSE`)

- returnInfo:

  Return all the data used in DFA (default = `FALSE`)

- silent:

  Silent-ish mode

- noTitle:

  Do not generate a title (only the subtitle)

- tsName:

  Name of y added as a subtitle to the plot

## Value

A dataframe with the Allan Factor (variance), Alan standard deviation
and error due to bin size

## See also

Other Fluctuation Analyses: [`fd_RR()`](fd_RR.md),
[`fd_dfa()`](fd_dfa.md), [`fd_mfdfa()`](fd_mfdfa.md),
[`fd_psd()`](fd_psd.md), [`fd_sda()`](fd_sda.md),
[`fd_sev()`](fd_sev.md)
