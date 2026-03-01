# Get (C)RQA measures from a Recurrence Matrix

Use `rp_measures`

## Usage

``` r
rp_measures_main(
  RM,
  emRad = NULL,
  DLmin = 2,
  VLmin = 2,
  HLmin = 2,
  DLmax = length(Matrix::diag(RM)),
  VLmax = length(Matrix::diag(RM)),
  HLmax = length(Matrix::diag(RM)),
  AUTO = NULL,
  theiler = NULL,
  chromatic = FALSE,
  matrices = FALSE,
  doHalf = FALSE,
  Nboot = NULL,
  CL = 0.95,
  doParallel = FALSE
)
```

## Arguments

- RM:

  A distance matrix (set `emRad = NA` to estimate a radius), or a matrix
  of zeroes and ones

- emRad:

  Threshold for distance value that counts as a recurrence (ignored is
  `RM` is a binary matrix)

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

- AUTO:

  Auto-recurrence? (default = `FALSE`)

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

- chromatic:

  Force chromatic RQA? If `NA` the value of the `RM` attribute
  `"chromatic"` will be used, if present (default = `NA`)

- matrices:

  Return matrices? (default = `FALSE`)

- Nboot:

  How many bootstrap replications? (default = `NULL`)

- CL:

  Confidence limit for bootstrap results (default = `.95`)

## See also

Other Recurrence Quantification Analysis: [`rp_cl()`](rp_cl.md),
[`rp_measures()`](rp_measures.md)
