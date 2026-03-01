# Line length distributions

Extract lengths of diagonal, vertical and horizontal line segments from
a recurrence matrix.

## Usage

``` r
rp_lineDist(
  RM,
  DLmin = 2,
  VLmin = 2,
  HLmin = 2,
  DLmax = length(Matrix::diag(RM)) - 1,
  VLmax = length(Matrix::diag(RM)) - 1,
  HLmax = length(Matrix::diag(RM)) - 1,
  d = NULL,
  theiler = NA,
  recurrenceTimes = FALSE,
  AUTO = NULL,
  chromatic = FALSE,
  matrices = FALSE
)
```

## Arguments

- RM:

  A distance matrix (set `emRad = NA` to estimate a radius), or a matrix
  of zeroes and ones

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

- d:

  Vector of diagonals to be extracted from matrix `RP` before line
  length distributions are calculated. A one element vector will be
  interpreted as a windowsize, e.g., `d = 50` will extract the diagonal
  band `-50:50`. A two element vector will be interpreted as a band,
  e.g. `d = c(-50,100)` will extract diagonals `-50:100`. If
  `length(d) > 2`, the numbers will be interpreted to refer to
  individual diagonals, `d = c(-50,50,100)` will extract diagonals
  `-50,50,100`. If `length(d)` is `NULL`, 1 or 2, the theiler window is
  applied before diagonals are extracted. The theiler window is ignored
  if `length(d)>2`, or if it is larger than the matrix or band indicated
  by parameter `d`. A warning will be given is a theiler window was
  already applied to the matrix.

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

- recurrenceTimes:

  Relevant for Recurrence Time analysis: Return the distribution of 0
  valued segments in nonzero diagonals/verticals/horizontals. This
  indicates the time between subsequent line structures.

- AUTO:

  Auto-recurrence? (default = `FALSE`)

- chromatic:

  Perform a chromatic RQA. This assumes the recurring values represent
  the labels of an unordered categorical variable (default = `FALSE`)

- matrices:

  Return matrices? (default = `FALSE`)

## Value

A list object with distributions of line lengths. If `matrices = TRUE`
datafr are returned whose columns represent the nonzero diagonals,
verticals, or, horizontals.

## Details

Based on the Matlab function `linedists` by Stefan Schinkel, Copyright
(C) 2009 Stefan Schinkel, University of Potsdam,
http://www.agnld.uni-potsdam.de

References: S. Schinkel, N. Marwan, O. Dimigen & J. Kurths (2009):
"Confidence Bounds of recurrence-based complexity measures Physics
Letters A, 373(26), pp. 2245-2250

Copyright (C) 2009 Stefan Schinkel, University of Potsdam
<http://www.agnld.uni-potsdam.de>

## See also

Other Distance matrix operations (recurrence plot):
[`bandReplace()`](bandReplace.md),
[`createCorridor()`](createCorridor.md), [`mat_di2bi()`](mat_di2bi.md),
[`mat_di2ch()`](mat_di2ch.md), [`mat_di2we()`](mat_di2we.md),
[`mat_hamming()`](mat_hamming.md), [`rp()`](rp.md),
[`rp_nzdiags()`](rp_nzdiags.md), [`rp_plot()`](rp_plot.md),
[`rp_size()`](rp_size.md)

## Author

Fred Hasselman
