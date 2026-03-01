# Corridor analysis

Create a corridor around the main diagonal. For long time series it may
not make sense to evaluate recurrences on the longest time scales.

## Usage

``` r
createCorridor(mat, lower, upper, value = NA, silent = TRUE)
```

## Arguments

- mat:

  A Matrix

- lower:

  Lower diagonal to be included in the corridor (should be \\\le 0\\)

- upper:

  Upper diagonal to be included in the corridor (should be \\\ge 0\\)

- value:

  A single value to replace all values outside the corridor (default =
  `NA`)

- silent:

  Operate in silence, only (some) warnings will be shown (default =
  `TRUE`)

## Value

A matrix in which the values outside the corridor have been replaced

## See also

Other Distance matrix operations (recurrence plot):
[`bandReplace()`](bandReplace.md), [`mat_di2bi()`](mat_di2bi.md),
[`mat_di2ch()`](mat_di2ch.md), [`mat_di2we()`](mat_di2we.md),
[`mat_hamming()`](mat_hamming.md), [`rp()`](rp.md),
[`rp_lineDist()`](rp_lineDist.md), [`rp_nzdiags()`](rp_nzdiags.md),
[`rp_plot()`](rp_plot.md), [`rp_size()`](rp_size.md)

## Author

Fred Hasselman

## Examples

``` r
# Create a 10 by 10 matrix
library(Matrix)
m <- Matrix(rnorm(10),10,10)

createCorridor(m,-7,7,0)   # Set diagonals 10 9 and 8 to 0.
#> 10 x 10 sparse Matrix of class "dgCMatrix"
#>                                                                              
#>  [1,] -0.13399701 -0.13399701 -0.13399701 -0.13399701 -0.13399701 -0.13399701
#>  [2,] -1.91008747 -1.91008747 -1.91008747 -1.91008747 -1.91008747 -1.91008747
#>  [3,] -0.27923724 -0.27923724 -0.27923724 -0.27923724 -0.27923724 -0.27923724
#>  [4,] -0.31344598 -0.31344598 -0.31344598 -0.31344598 -0.31344598 -0.31344598
#>  [5,]  1.06730788  1.06730788  1.06730788  1.06730788  1.06730788  1.06730788
#>  [6,]  0.07003485  0.07003485  0.07003485  0.07003485  0.07003485  0.07003485
#>  [7,] -0.63912332 -0.63912332 -0.63912332 -0.63912332 -0.63912332 -0.63912332
#>  [8,] -0.04996490 -0.04996490 -0.04996490 -0.04996490 -0.04996490 -0.04996490
#>  [9,] -0.25148344 -0.25148344 -0.25148344 -0.25148344 -0.25148344 -0.25148344
#> [10,]  0.44479712  0.44479712  0.44479712  0.44479712  0.44479712  0.44479712
#>                                                      
#>  [1,] -0.13399701 -0.13399701 -0.13399701 -0.13399701
#>  [2,] -1.91008747 -1.91008747 -1.91008747 -1.91008747
#>  [3,] -0.27923724 -0.27923724 -0.27923724 -0.27923724
#>  [4,] -0.31344598 -0.31344598 -0.31344598 -0.31344598
#>  [5,]  1.06730788  1.06730788  1.06730788  1.06730788
#>  [6,]  0.07003485  0.07003485  0.07003485  0.07003485
#>  [7,] -0.63912332 -0.63912332 -0.63912332 -0.63912332
#>  [8,] -0.04996490 -0.04996490 -0.04996490 -0.04996490
#>  [9,] -0.25148344 -0.25148344 -0.25148344 -0.25148344
#> [10,]  0.44479712  0.44479712  0.44479712  0.44479712
```
