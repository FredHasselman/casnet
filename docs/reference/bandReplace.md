# Replace matrix diagonals

Sets a band of matrix diagonals to any given value

## Usage

``` r
bandReplace(mat, lower, upper, value = NA, silent = TRUE)
```

## Arguments

- mat:

  A Matrix

- lower:

  Lower diagonal to be included in the band (should be \\\le 0\\)

- upper:

  Upper diagonal to be included in the band (should be \\\ge 0\\)

- value:

  A single value to replace all values in the selected band (default =
  `NA`)

- silent:

  Operate in silence, only (some) warnings will be shown (default =
  `TRUE`)

## Value

A matrix in which the values in the selected diagonals have been
replaced

## See also

Other Distance matrix operations (recurrence plot):
[`createCorridor()`](createCorridor.md), [`mat_di2bi()`](mat_di2bi.md),
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

bandReplace(m,-1,1,0)   # Replace diagonal and adjacent bands with 0 (Theiler window of 1)
#> 'as(<dgeMatrix>, "dgCMatrix")' is deprecated.
#> Use 'as(., "CsparseMatrix")' instead.
#> See help("Deprecated") and help("Matrix-deprecated").
#> 10 x 10 sparse Matrix of class "dgCMatrix"
#>                                                                           
#>  [1,]  .           .           1.8885049  1.8885049  1.8885049  1.88850493
#>  [2,]  .           .           .         -0.0974451 -0.0974451 -0.09744510
#>  [3,] -0.93584735  .           .          .         -0.9358474 -0.93584735
#>  [4,] -0.01595031 -0.01595031  .          .          .         -0.01595031
#>  [5,] -0.82678895 -0.82678895 -0.8267890  .          .          .         
#>  [6,] -1.51239965 -1.51239965 -1.5123997 -1.5123997  .          .         
#>  [7,]  0.93536319  0.93536319  0.9353632  0.9353632  0.9353632  .         
#>  [8,]  0.17648861  0.17648861  0.1764886  0.1764886  0.1764886  0.17648861
#>  [9,]  0.24368546  0.24368546  0.2436855  0.2436855  0.2436855  0.24368546
#> [10,]  1.62354888  1.62354888  1.6235489  1.6235489  1.6235489  1.62354888
#>                                                      
#>  [1,]  1.88850493  1.88850493  1.88850493  1.88850493
#>  [2,] -0.09744510 -0.09744510 -0.09744510 -0.09744510
#>  [3,] -0.93584735 -0.93584735 -0.93584735 -0.93584735
#>  [4,] -0.01595031 -0.01595031 -0.01595031 -0.01595031
#>  [5,] -0.82678895 -0.82678895 -0.82678895 -0.82678895
#>  [6,]  .          -1.51239965 -1.51239965 -1.51239965
#>  [7,]  .           .           0.93536319  0.93536319
#>  [8,]  .           .           .           0.17648861
#>  [9,]  0.24368546  .           .           .         
#> [10,]  1.62354888  1.62354888  .           .         
```
