# Distance to binary matrix

Distance matrix to binary matrix based on threshold value

## Usage

``` r
mat_di2bi(distmat, emRad = NA, convMat = FALSE)
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

A (sparse) matrix with only 0s and 1s

## See also

Other Distance matrix operations (recurrence plot):
[`bandReplace()`](bandReplace.md),
[`createCorridor()`](createCorridor.md), [`mat_di2ch()`](mat_di2ch.md),
[`mat_di2we()`](mat_di2we.md), [`mat_hamming()`](mat_hamming.md),
[`rp()`](rp.md), [`rp_lineDist()`](rp_lineDist.md),
[`rp_nzdiags()`](rp_nzdiags.md), [`rp_plot()`](rp_plot.md),
[`rp_size()`](rp_size.md)

Other Distance matrix operations (recurrence network):
[`mat_di2ch()`](mat_di2ch.md), [`mat_di2we()`](mat_di2we.md),
[`rn()`](rn.md), [`rn_phaseInfo()`](rn_phaseInfo.md),
[`rn_phases()`](rn_phases.md), [`rn_plot()`](rn_plot.md),
[`rn_recSpec()`](rn_recSpec.md)
