# Weighted to Binary matrix

Weighted to Binary matrix

## Usage

``` r
mat_we2bi(distmat, emRad, theiler = 0, convMat = FALSE)
```

## Arguments

- distmat:

  Distance matrix

- emRad:

  The radius or threshold value

- convMat:

  Should the matrix be converted from a `distmat` object of class
  [`Matrix::Matrix()`](https://rdrr.io/pkg/Matrix/man/Matrix.html) to
  [`base::matrix()`](https://rdrr.io/r/base/matrix.html) (or vice versa)

## Value

A binary matrix
