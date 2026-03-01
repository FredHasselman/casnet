# Cross Diagonal Recurrence Profiles

Cross Diagonal Recurrence Profiles

## Usage

``` r
rqa_diagProfile2(
  y1,
  y2 = NULL,
  emDim = 1,
  emLag = 1,
  emRad = NULL,
  targetValue = 0.05,
  theiler = 0,
  corridor = NA,
  doEmbed = TRUE,
  chromatic = FALSE,
  AUTO = FALSE,
  method = c("Euclidean", "SBD")[1],
  doPlot = TRUE,
  silent = TRUE
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

- targetValue:

  A value passed to `est_radius(...,type="fixed", targetMeasure="RR")`
  if `is.na(emRad)==TRUE`.

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

- doEmbed:

  If `FALSE`, a distance matrix will be returned that is not embedded by
  `emDim` and `emLag` (Multidimensional RQA). If `y1` and/or `y2` are
  data frames, the columns will be used as the state space dimensions
  (default = `TRUE`)

- chromatic:

  Perform a chromatic RQA. This assumes the recurring values represent
  the labels of an unordered categorical variable (default = `FALSE`)

- AUTO:

  Auto-recurrence? (default = `FALSE`)

- method:

  Distance measure to use. Any option that is valid for argument
  `method` of `proxy::dist()`. Type `proxy::pr_DB$get_entries()` to see
  a list of all the options. Common methods are:
  `"Euclidean", "Manhattan", "Minkowski", "Chebysev"` (or the same but
  shorter: `"L2","L1","Lp", "max"` distance). To use the shape based
  distance for phase-based recurrence use `"SBD"` (default =
  `"Euclidean"`)

- doPlot:

  Plot the matrix by calling [`rp_plot()`](rp_plot.md) with default
  settings

- silent:

  Silent-ish mode

## Value

Data file with CDRP output
