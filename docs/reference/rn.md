# Create a Recurrence Network Matrix

This function serves as a wrapper for function [`rp()`](rp.md), it will
add some attributes to the matrix related to network representation.
These attributes will be used to decide which network type to generate
(e.g. undirected, directed, weighted, etc.)

## Usage

``` r
rn(
  y1,
  y2 = NULL,
  emDim = 1,
  emLag = 1,
  emRad = NULL,
  theiler = 0,
  directed = FALSE,
  cumulative = TRUE,
  weighted = FALSE,
  weightedBy = c("none", "si", "rt", "rf")[1],
  rescaleWeights = FALSE,
  fs = NA,
  to.ts = NULL,
  order.by = NULL,
  to.sparse = FALSE,
  method = c("Euclidean", "max", "SBD")[1],
  rescaleDist = c("none", "maxDist", "meanDist")[1],
  targetValue = 0.05,
  returnGraph = FALSE,
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

- directed:

  Should the matrix be considered to represent a directed network?
  (default = `FALSE`)

- cumulative:

  To make the network represent cumulative time, set `directed = TRUE`
  and `cumulative = TRUE`. This will set the upper triangle of the
  recurrence matrix to `0` and ensures that the network edges represent
  recurrent values that have occurred in the `past` relative to the
  current observed value (node). If `directed = FALSE` the argument is
  ignored (default = `TRUE`).

- weighted:

  Should the matrix be considered to represent a weighted network?
  (default = `FALSE`)

- weightedBy:

  After setting values smaller than `emRad` to `0`, what should the
  recurrent values represent? The default is to use the state space
  similarity (distance/proximity) values as weights (`"si"`). Other
  option are `"rt"` for *recurrence time* and `"rf"` for *recurrence
  time frequency*, Because vertices represent time points in
  \\\epsilon\\-thresholded recurrence networks, a difference of two
  vertex-indices represents duration. If an edge `e1` connects `v1` and
  `v10` then the *recurrence time* will be the difference of the vertex
  indices, `9`, and the *recurrence time frequency* will be `1/9`.

- rescaleWeights:

  If set to `TRUE` and `weighted = TRUE`, all weight values will be
  rescaled to `[0,1]`, where `0` means no recurrence relation and `1`
  the maximum weight value.

- fs:

  Sample frequency: A numeric value interpreted as the
  `number of observed samples per unit of time`. If the weights
  represent recurrence times (`"rt"`), they will be divided by the value
  in `fs`. If the weights represent recurrence time frequencies
  (`"rf"`), they will be multiplied by the value of `fs` (default =
  `NA`)

- to.ts:

  Should `y1` and `y2` be converted to time series objects?

- order.by:

  If `to.ts = TRUE`, pass a vector of the same length as `y1` and `y2`.
  It will be used as the time index, if `NA` the vector indices will be
  used to represent time.

- to.sparse:

  Should sparse matrices be used?

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

- returnGraph:

  Return an
  [`igraph::igraph()`](https://r.igraph.org/reference/aaa-igraph-package.html)
  object (default = `FALSE`)

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

  Any paramters to pass to [`rn_plot()`](rn_plot.md) if `doPlot = TRUE`

## Value

A (Coss-) Recurrence matrix that can be interpreted as an adjacency (or
incidence) matrix.

## See also

Other Distance matrix operations (recurrence network):
[`mat_di2bi()`](mat_di2bi.md), [`mat_di2ch()`](mat_di2ch.md),
[`mat_di2we()`](mat_di2we.md), [`rn_phaseInfo()`](rn_phaseInfo.md),
[`rn_phases()`](rn_phases.md), [`rn_plot()`](rn_plot.md),
[`rn_recSpec()`](rn_recSpec.md)
