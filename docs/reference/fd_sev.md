# Calculate FD using Sevcik's method

Calculate FD using Sevcik's method

## Usage

``` r
fd_sev(
  y,
  detrend = FALSE,
  adjustSumOrder = FALSE,
  smallNapprox = FALSE,
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

  A time series or numeric vector

- detrend:

  Subtract linear trend from the series (default = `TRUE`).

- adjustSumOrder:

  Adjust the time series (summation or differencing), based on the
  global scaling exponent, see e.g.
  <https://www.frontiersin.org/files/Articles/23948/fphys-03-00141-r2/image_m/fphys-03-00141-t001.jpg>Ihlen (2012)
  (default = `TRUE`)

- smallNapprox:

  Force use of small sample approximation (default for N \< 128)

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

An FD estimate

## References

Sevcik, C. (1998). A procedure to Estimate the Fractal Dimension of
Waveforms. Paper available at http://arxiv.org/pdf/1003.5266.pdf

## See also

Other Fluctuation Analyses: [`fd_RR()`](fd_RR.md),
[`fd_allan()`](fd_allan.md), [`fd_dfa()`](fd_dfa.md),
[`fd_mfdfa()`](fd_mfdfa.md), [`fd_psd()`](fd_psd.md),
[`fd_sda()`](fd_sda.md)

## Author

Fred Hasselman
