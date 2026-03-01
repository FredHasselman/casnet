# rp_nzdiags

Get all nonzero diagonals of a binary matrix, or, diagonals specified as
a vector by argument `d`.

## Usage

``` r
rp_nzdiags(
  RM = NULL,
  d = NULL,
  returnVectorList = TRUE,
  returnNZtriplets = FALSE,
  removeNZ = TRUE,
  silent = TRUE
)
```

## Arguments

- RM:

  A binary (0,1) matrix.

- d:

  An optional vector of diagonals to extract.

- returnVectorList:

  Return list

- returnNZtriplets:

  Return a dataframe with coordinates of only nonzero elements in
  diagonals (default = `FALSE`)

- removeNZ:

  Remove nonzero diagonals if `TRUE`. If `FALSE` returns the full
  diagonals matrix. Use e.g. to plot diagonal recurrence profiles
  (default = `TRUE`)

- silent:

  Silent-ish mode

## Value

A matrix object with nonzero diagonals as columns and/or a dataframe
with coordinates of nonzero diagonal elements

## See also

Other Distance matrix operations (recurrence plot):
[`bandReplace()`](bandReplace.md),
[`createCorridor()`](createCorridor.md), [`mat_di2bi()`](mat_di2bi.md),
[`mat_di2ch()`](mat_di2ch.md), [`mat_di2we()`](mat_di2we.md),
[`mat_hamming()`](mat_hamming.md), [`rp()`](rp.md),
[`rp_lineDist()`](rp_lineDist.md), [`rp_plot()`](rp_plot.md),
[`rp_size()`](rp_size.md)

## Author

Fred Hasselman
