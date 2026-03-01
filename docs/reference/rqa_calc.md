# rqa_calc

rqa_calc

## Usage

``` r
rqa_calc(
  y1,
  y2 = NA,
  emRad = NULL,
  DLmin = 2,
  VLmin = 2,
  HLmin = 2,
  DLmax = NA,
  VLmax = NA,
  HLmax = NA,
  theiler = NA,
  corridor = NA,
  AUTO = NULL,
  method = "Euclidean",
  includeDiagonal = NA,
  chromatic = FALSE,
  anisotropyHV = FALSE,
  asymmetryUL = FALSE,
  returnUL = FALSE,
  recurrenceTimes = FALSE,
  distributions = FALSE
)
```

## Arguments

- y1:

  A numeric vector or time series

- y2:

  A numeric vector or time series for cross recurrence

- emRad:

  The threshold (emRad) to apply to the distance matrix to create a
  binary or weighted matrix. If `NULL`, an unthresholded matrix will be
  created (default = `NULL`)

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

- corridor:

  corridor

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

- chromatic:

  Perform a chromatic RQA. This assumes the recurring values represent
  the labels of an unordered categorical variable (default = `FALSE`)

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

- returnUL:

  Return the (C)RQA values for the upper and lower triangle on which
  asymmetry ratio calculations are based? (default = `FALSE`)

- recurrenceTimes:

  Return measures based on 'white lines', the recurrence times (default
  = `FALSE`)

- distributions:

  Return line length distributions.

## Value

CRQA measures and matrices of line distributions (if requested) based on
massively parallel analysis, without building the recurrence matrix.
