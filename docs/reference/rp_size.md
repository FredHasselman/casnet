# rp_size

Calculate the maximum possible number of recurrent points in a
recurrence matrix.

## Usage

``` r
rp_size(RM = NULL, dims = NULL, AUTO = NULL, theiler = NULL)
```

## Arguments

- RM:

  A Matrix object

- dims:

  Two element vector representing the dimensions of Matrix `RM`. If
  `dims` is provided, the Matrix does not have to be passed as an
  argument (default = `NA`)

- AUTO:

  Is the Matrix an Auto Recurrence Matrix? If so, the length of the
  diagonal will be subtracted from the matrix size, pass `FALSE` to
  prevent this behaviour. If `NULL` (default) `AUTO` will take on the
  value of `isSymmetric(RM)`.

- theiler:

  Should a Theiler window be applied?

## Value

Matrix size for computation of recurrence measures.

## Details

This function can take into account the presence of a `theiler` window,
that is the points in the window will be excluded from the calculation.
For example, some scholars will exclude the main diagonal from the
calculation of the recurrence rate.

## See also

Other Distance matrix operations (recurrence plot):
[`bandReplace()`](bandReplace.md),
[`createCorridor()`](createCorridor.md), [`mat_di2bi()`](mat_di2bi.md),
[`mat_di2ch()`](mat_di2ch.md), [`mat_di2we()`](mat_di2we.md),
[`mat_hamming()`](mat_hamming.md), [`rp()`](rp.md),
[`rp_lineDist()`](rp_lineDist.md), [`rp_nzdiags()`](rp_nzdiags.md),
[`rp_plot()`](rp_plot.md)

## Examples

``` r
# Create a 10 by 10 matrix
library(Matrix)
m <- Matrix(rnorm(10),10,10)

rp_size(RM = m, AUTO = TRUE, theiler = 0)  # Subtract diagonal
#> $rp_size_total
#> [1] 100
#> 
#> $rp_size_theiler
#> [1] 100
#> 
rp_size(RM = m, AUTO = FALSE,theiler = 0)  # Do not subtract diagonal
#> $rp_size_total
#> [1] 100
#> 
#> $rp_size_theiler
#> [1] 100
#> 
rp_size(RM = m, AUTO = NULL, theiler = 0)  # Matrix is symmetrical, AUTO is set to TRUE
#> $rp_size_total
#> [1] 100
#> 
#> $rp_size_theiler
#> [1] 100
#> 
rp_size(RM = m, AUTO = NULL, theiler = 1)  # Subtract a Theiler window of 1 around and including the diagonal
#> $rp_size_total
#> [1] 100
#> 
#> $rp_size_theiler
#> [1] 90
#> 

# Calculate without a matrix
rp_size(dims = c(10,10), AUTO = TRUE, TRUE,0)
#> $rp_size_total
#> [1] 100
#> 
#> $rp_size_theiler
#> [1] 100
#> 
rp_size(dims = c(10,10), AUTO = FALSE,FALSE,0)
#> $rp_size_total
#> [1] 100
#> 
#> $rp_size_theiler
#> [1] 100
#> 
```
