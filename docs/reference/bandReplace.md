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
#>  [1,]  .           .          -0.91407483 -0.91407483 -0.91407483 -0.91407483
#>  [2,]  .           .           .           0.46815442  0.46815442  0.46815442
#>  [3,]  0.36295126  .           .           .           0.36295126  0.36295126
#>  [4,] -1.30454355 -1.30454355  .           .           .          -1.30454355
#>  [5,]  0.73777632  0.73777632  0.73777632  .           .           .         
#>  [6,]  1.88850493  1.88850493  1.88850493  1.88850493  .           .         
#>  [7,] -0.09744510 -0.09744510 -0.09744510 -0.09744510 -0.09744510  .         
#>  [8,] -0.93584735 -0.93584735 -0.93584735 -0.93584735 -0.93584735 -0.93584735
#>  [9,] -0.01595031 -0.01595031 -0.01595031 -0.01595031 -0.01595031 -0.01595031
#> [10,] -0.82678895 -0.82678895 -0.82678895 -0.82678895 -0.82678895 -0.82678895
#>                                                   
#>  [1,] -0.91407483 -0.9140748 -0.9140748 -0.9140748
#>  [2,]  0.46815442  0.4681544  0.4681544  0.4681544
#>  [3,]  0.36295126  0.3629513  0.3629513  0.3629513
#>  [4,] -1.30454355 -1.3045435 -1.3045435 -1.3045435
#>  [5,]  0.73777632  0.7377763  0.7377763  0.7377763
#>  [6,]  .           1.8885049  1.8885049  1.8885049
#>  [7,]  .           .         -0.0974451 -0.0974451
#>  [8,]  .           .          .         -0.9358474
#>  [9,] -0.01595031  .          .          .        
#> [10,] -0.82678895 -0.8267890  .          .        
```
