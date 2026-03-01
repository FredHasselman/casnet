# Plot (thresholded) distance matrix as a network

Plot (thresholded) distance matrix as a network

## Usage

``` r
rn_plot(
  RN,
  plotDimensions = FALSE,
  plotMeasures = FALSE,
  drawGrid = FALSE,
  markEpochsLOI = NULL,
  radiusValue = NA,
  title = "",
  xlabel = "",
  ylabel = "",
  plotSurrogate = NA,
  returnOnlyObject = FALSE
)
```

## Arguments

- RN:

  A distance matrix or recurrence matrix

- plotDimensions:

  Should the state vectors be plotted if they are available as
  attributes of RM (default = `TRUE`)

- plotMeasures:

  Print common (C)RQA measures in the plot if the matrix is binary
  (default = `FALSE`)

- drawGrid:

  Draw a grid on the recurrence plot (default = `FALSE`)

- markEpochsLOI:

  Pass a factor whose levels indicate different epochs or phases in the
  time series and use the line of identity to represent the levels by
  different colours (default = `NULL`)

- radiusValue:

  If `plotMeasures = TRUE` and RM is an unthresholded matrix, this value
  will be used to calculate recurrence measures. If
  `plotMeasures = TRUE` and RM is already a binary recurrence matrix,
  pass the radius that was used as a threshold to create the matrix for
  display purposes. If `plotMeasures = TRUE` and `radiusValue = NA`,
  function [`est_radius()`](est_radius.md) will be called with default
  settings (find a radius that yields `.05` recurrence rate). If
  `plotMeasures = FALSE` this setting will be ignored.

- title:

  A title for the plot

- xlabel:

  An x-axis label

- ylabel:

  An y-axis label

- plotSurrogate:

  Should a 2-panel comparison plot based on surrogate time series be
  added? If `RM` has attributes `y1` and `y2` containing the time series
  data (i.e. it was created by a call to [rp](rp.md)), the following
  options are available: "RS" (random shuffle), "RP" (randomised
  phases), "AAFT" (amplitude adjusted fourier transform). If no
  timeseries data is included, the columns will be shuffled. NOTE: This
  is not a surrogate test, just 1 surrogate is created from `y1`.
  (default = `FALSE`)

- returnOnlyObject:

  Return the ggplot object only, do not draw the plot (default = `TRUE`)

## Value

A nice plot of the recurrence network

## See also

Other Distance matrix operations (recurrence network):
[`mat_di2bi()`](mat_di2bi.md), [`mat_di2ch()`](mat_di2ch.md),
[`mat_di2we()`](mat_di2we.md), [`rn()`](rn.md),
[`rn_phaseInfo()`](rn_phaseInfo.md), [`rn_phases()`](rn_phases.md),
[`rn_recSpec()`](rn_recSpec.md)
