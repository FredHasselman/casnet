# Create a Distance Matrix

Create a Distance Matrix

## Usage

``` r
rp(
  y1,
  y2 = NULL,
  emDim = 1,
  emLag = 1,
  emRad = NULL,
  theiler = NA,
  to.ts = NULL,
  order.by = NULL,
  to.sparse = TRUE,
  weighted = FALSE,
  weightedBy = "si",
  method = c("Euclidean", "max", "SBD")[1],
  rescaleDist = c("none", "maxDist", "meanDist")[1],
  targetValue = 0.05,
  chromatic = FALSE,
  returnMeasures = FALSE,
  doPlot = FALSE,
  doEmbed = TRUE,
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

  The threshold (emRad) to apply to the distance matrix to create a
  binary or weighted matrix. If `NULL`, an unthresholded matrix will be
  created (default = `NULL`)

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

- to.ts:

  Should `y1` and `y2` be converted to time series objects?

- order.by:

  If `to.ts = TRUE`, pass a vector of the same length as `y1` and `y2`.
  It will be used as the time index, if `NA` the vector indices will be
  used to represent time.

- to.sparse:

  Should sparse matrices be used?

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

- returnMeasures:

  Should the output of [`rp_measures()`](rp_measures.md) be returned as
  an attribute `"measures"` to the matrix? If `silent = FALSE` results
  will also be output to the console. (default = `FALSE`)

- doPlot:

  Plot the matrix by calling [`rp_plot()`](rp_plot.md) with default
  settings

- doEmbed:

  If `FALSE`, a distance matrix will be returned that is not embedded by
  `emDim` and `emLag` (Multidimensional RQA). If `y1` and/or `y2` are
  data frames, the columns will be used as the state space dimensions
  (default = `TRUE`)

- silent:

  Silent-ish mode

- ...:

  Any parameters to pass to [`rp_plot()`](rp_plot.md) if `doPlot = TRUE`

## Value

A (Coss-) Recurrence matrix with attributes:

- `emdims1` and `emdims2` - A matrix of surrogate dimensions

- `emdims1.name` and `emdims2.name` - Names of surrogate dimensions

- `method` and `call` - The distance `method` used by `proxy::dist()`

- `weighted` - Whether a weighted matrix is returned

- `emDim`, `emLag` and `emRad` - The embedding parameters

- `AUTO` - Whether the matrix represents AUTO recurrence

## Note

The calculation of the (C)RQA measures in [casnet](casnet-package.md)
can be different from other packages. For example, depending on the
value of `theiler` the main diagonal can be included or excluded from
the calculations, whereas some software will always include the
diagonal.

## See also

Other Distance matrix operations (recurrence plot):
[`bandReplace()`](bandReplace.md),
[`createCorridor()`](createCorridor.md), [`mat_di2bi()`](mat_di2bi.md),
[`mat_di2ch()`](mat_di2ch.md), [`mat_di2we()`](mat_di2we.md),
[`mat_hamming()`](mat_hamming.md), [`rp_lineDist()`](rp_lineDist.md),
[`rp_nzdiags()`](rp_nzdiags.md), [`rp_plot()`](rp_plot.md),
[`rp_size()`](rp_size.md)
