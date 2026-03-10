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
#>  [1,] -1.5123997 -1.5123997 -1.5123997 -1.5123997 -1.5123997 -1.5123997
#>  [2,]  0.9353632  0.9353632  0.9353632  0.9353632  0.9353632  0.9353632
#>  [3,]  0.1764886  0.1764886  0.1764886  0.1764886  0.1764886  0.1764886
#>  [4,]  0.2436855  0.2436855  0.2436855  0.2436855  0.2436855  0.2436855
#>  [5,]  1.6235489  1.6235489  1.6235489  1.6235489  1.6235489  1.6235489
#>  [6,]  0.1120381  0.1120381  0.1120381  0.1120381  0.1120381  0.1120381
#>  [7,] -0.1339970 -0.1339970 -0.1339970 -0.1339970 -0.1339970 -0.1339970
#>  [8,] -1.9100875 -1.9100875 -1.9100875 -1.9100875 -1.9100875 -1.9100875
#>  [9,] -0.2792372 -0.2792372 -0.2792372 -0.2792372 -0.2792372 -0.2792372
#> [10,] -0.3134460 -0.3134460 -0.3134460 -0.3134460 -0.3134460 -0.3134460
#>                                                  
#>  [1,] -1.5123997 -1.5123997 -1.5123997 -1.5123997
#>  [2,]  0.9353632  0.9353632  0.9353632  0.9353632
#>  [3,]  0.1764886  0.1764886  0.1764886  0.1764886
#>  [4,]  0.2436855  0.2436855  0.2436855  0.2436855
#>  [5,]  1.6235489  1.6235489  1.6235489  1.6235489
#>  [6,]  0.1120381  0.1120381  0.1120381  0.1120381
#>  [7,] -0.1339970 -0.1339970 -0.1339970 -0.1339970
#>  [8,] -1.9100875 -1.9100875 -1.9100875 -1.9100875
#>  [9,] -0.2792372 -0.2792372 -0.2792372 -0.2792372
#> [10,] -0.3134460 -0.3134460 -0.3134460 -0.3134460
```
