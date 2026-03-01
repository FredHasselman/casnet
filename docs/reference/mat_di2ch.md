# Distance to chromatic matrix

Distance matrix to chromatic matrix based on unordered categorical
series

## Usage

``` r
mat_di2ch(distmat, y, emRad, theiler = 0, convMat = FALSE)
```

## Arguments

- distmat:

  Distance matrix

- y:

  One of the dimensions (as a data frame or matrix) of the RP which must
  contain unique unordered categorical values

- emRad:

  The radius or threshold value

- theiler:

  Use a theiler window around the line of identity / synchronisation to
  remove high auto-correlation at short time-lags (default = `0`)

- convMat:

  convMat Should the matrix be converted from a `distmat` object of
  class [`Matrix::Matrix()`](https://rdrr.io/pkg/Matrix/man/Matrix.html)
  to [`base::matrix()`](https://rdrr.io/r/base/matrix.html) (or vice
  versa)

## Value

A matrix with 0s and the unordered categorical values that are recurring

## See also

Other Distance matrix operations (recurrence plot):
[`bandReplace()`](bandReplace.md),
[`createCorridor()`](createCorridor.md), [`mat_di2bi()`](mat_di2bi.md),
[`mat_di2we()`](mat_di2we.md), [`mat_hamming()`](mat_hamming.md),
[`rp()`](rp.md), [`rp_lineDist()`](rp_lineDist.md),
[`rp_nzdiags()`](rp_nzdiags.md), [`rp_plot()`](rp_plot.md),
[`rp_size()`](rp_size.md)

Other Distance matrix operations (recurrence network):
[`mat_di2bi()`](mat_di2bi.md), [`mat_di2we()`](mat_di2we.md),
[`rn()`](rn.md), [`rn_phaseInfo()`](rn_phaseInfo.md),
[`rn_phases()`](rn_phases.md), [`rn_plot()`](rn_plot.md),
[`rn_recSpec()`](rn_recSpec.md)
