# Massively Parallel RQA analysis

**\[experimental\]** Calculate (C)RQA measures without creating a
recurrence matrix. Can handle very large time series and requires
package
[future.apply::future.apply](https://future.apply.futureverse.org/reference/future.apply.html)
to be installed.

## Usage

``` r
rqa_par(
  y1,
  y2 = NULL,
  emDim = 1,
  emLag = 1,
  emRad = NULL,
  theiler = NA,
  includeDiagonal = NA,
  AUTO = NULL,
  DLmin = 2,
  VLmin = 2,
  HLmin = 2,
  DLmax = NA,
  VLmax = NA,
  HLmax = NA,
  weighted = FALSE,
  weightedBy = "si",
  method = c("Euclidean", "SBD")[1],
  rescaleDist = c("none", "maxDist", "meanDist")[1],
  targetValue = 0.05,
  chromatic = FALSE,
  recurrenceTimes = FALSE,
  doPlot = FALSE,
  doEmbed = TRUE,
  anisotropyHV = FALSE,
  asymmetryUL = FALSE,
  silent = TRUE,
  ...
)
```

## Arguments

- y1:

  A numeric vector or time series

- y2:

  A numeric vector or time series for cross recurrence

- emDim:

  The embedding dimensions

- emLag:

  The embedding lag

- emRad:

  The radius threshold on the distance measure that determines the
  global Recurrence Rate. If the radius is not provided it will be
  estimated with default settings using function
  [est_radius](est_radius.md) (default = `NULL`)

- theiler:

  Use a `theiler` window around the main diagonal (Line of
  Identity/Synchronisation) to remove auto-correlations at short
  time-lags:

  - `0` will include the main diagonal in all RQA measure calculations.

  - `1` will remove the main diagonal from all RQA measure calculations.

  - `NA` (default), will check if the matrix is symmetrical , if so, it
    will remove the diagonal by setting `theiler = 1` (Line of Identity,
    Auto-RQA), if it is not symmetrical (Line of Synchronisation,
    Cross-RQA) it will set `theiler = 0`.

  - A value greater than `1` will remove that many diagonals around and
    including the diagonal from all RQA measure calculations. So
    `theiler = 2` means exclude `2` diagonals around the main diagonal,
    including the main diagonal itself: `[-1,0,1]`. If `theiler` is a
    numeric vector of `length(theiler) == 2` it is possible to exclude
    an asymmetrical window. The values are interpreted as end points in
    a sequence of diagonal ID's, e.g. `theiler = c(-1,5)` will exclude
    `[-1,0,1,2,3,4,5]`. If `length(theiler) > 2`, the values will be
    considered individual diagonal ID's, e.g.
    `theiler = c(-3,-1,0,2,5)`, will exclude only those specific ID's.
    Also see the note.

- AUTO:

  Auto-recurrence? (default = `FALSE`)

- DLmin:

  Minimal diagonal line length (default = `2`)

- VLmin:

  Minimal vertical line length (default = `2`)

- HLmin:

  Minimal horizontal line length (default = `2`)

- DLmax:

  Maximal diagonal line length (default = length of diagonal -1)

- VLmax:

  Maximal vertical line length (default = length of diagonal -1)

- HLmax:

  Maximal horizontal line length (default = length of diagonal -1)

- weighted:

  If `FALSE` a binary matrix will be returned. If `TRUE` every value
  larger than `emRad` will be `0`, but values smaller than `emRad` will
  be retained (default = `FALSE`)

- weightedBy:

  After setting values smaller than `emRad` to `0`, what should the
  recurrent values represent? The default is to use the state space
  similarity (distance/proximity) values as weights (`"si"`). Other
  option are `"rt"` for *recurrence time* and `"rf"` for *recurrence
  time frequency* (default = `"si"`)

- method:

  Distance measure to use. Any option that is valid for argument
  `method` of `proxy::dist()`. Type `proxy::pr_DB$get_entries()` to see
  a list of all the options. Common methods are:
  `"Euclidean", "Manhattan", "Minkowski", "Chebysev"` (or the same but
  shorter: `"L2","L1","Lp", "max"` distance). To use the shape based
  distance for phase-based recurrence use `"SBD"` (default =
  `"Euclidean"`)

- rescaleDist:

  Should the distance matrix be rescaled? Options are "none", "maxDist"
  to create a unit scale, "meanScale" to creat z-scores based on the
  mean distance. (default = `"none"`)

- targetValue:

  A value passed to `est_radius(...,type="fixed", targetMeasure="RR")`
  if `is.na(emRad)==TRUE`.

- chromatic:

  Perform a chromatic RQA. This assumes the recurring values represent
  the labels of an unordered categorical variable (default = `FALSE`)

- recurrenceTimes:

  Return measures based on 'white lines', the recurrence times (default
  = `FALSE`)

- doPlot:

  Plot the matrix by calling [`rp_plot()`](rp_plot.md) with default
  settings

- doEmbed:

  If `FALSE`, a distance matrix will be returned that is not embedded by
  `emDim` and `emLag` (Multidimensional RQA). If `y1` and/or `y2` are
  data frames, the columns will be used as the state space dimensions
  (default = `TRUE`)

- anisotropyHV:

  Return anisotropy ratio measures based on Horizontal and Vertical
  lines. The ratios are calculated as
  `(horizontal - vertical) / (horizontal + vertical)`. So a value of 0
  means no anisotropy, negative ratios indicate the measures based on
  vertical lines had higher values, positive ratios indicate the
  measures based on horizontal lines had higher values (default =
  `FALSE`)

- asymmetryUL:

  Return asymmetry ratio measures based on Upper and Lower triangles.
  The ratios are calculated as `(upper - lower) / (upper + lower)`. So a
  value of 0 means no asymmetry, negative ratios indicate the measures
  based on the lower triangle had the higher values, positive ratios
  indicate measures based on the upper triangle had higher values
  (default = `FALSE`)

- silent:

  Silent-ish mode

- ...:

  Any parameters to pass to [`rp_plot()`](rp_plot.md) if `doPlot = TRUE`

## Value

RQA
