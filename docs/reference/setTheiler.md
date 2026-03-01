# Set theiler window on a distance matrix or recurrence matrix.

Set theiler window on a distance matrix or recurrence matrix.

## Usage

``` r
setTheiler(RM, theiler = NA, silent = FALSE, chromatic = FALSE)
```

## Arguments

- RM:

  A distance matrix (set `emRad = NA` to estimate a radius), or a matrix
  of zeroes and ones

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

- silent:

  Silent-ish mode

- chromatic:

  Perform a chromatic RQA. This assumes the recurring values represent
  the labels of an unordered categorical variable (default = `FALSE`)

## Value

The matrix with the diagonals indicated in the `theiler` argument set to
either `max(RM)+1` (if `RM` is a distance matrix) or `0` (if `RM` is a
recurrence matrix).
