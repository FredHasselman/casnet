# Calculate Hamming distance

Calculate Hamming distance

## Usage

``` r
mat_hamming(X, Y = NULL, embedded = TRUE)
```

## Arguments

- X:

  A matrix (of coordinates)

- Y:

  A matrix (of coordinates)

- embedded:

  Do X and/or Y represent surrogate dimensions of an embedded time
  series?

## Value

A hamming-distance matrix of X, or X and Y. Useful for ordered and
unordered categorical data.

## See also

Other Distance matrix operations (recurrence plot):
[`bandReplace()`](bandReplace.md),
[`createCorridor()`](createCorridor.md), [`mat_di2bi()`](mat_di2bi.md),
[`mat_di2ch()`](mat_di2ch.md), [`mat_di2we()`](mat_di2we.md),
[`rp()`](rp.md), [`rp_lineDist()`](rp_lineDist.md),
[`rp_nzdiags()`](rp_nzdiags.md), [`rp_plot()`](rp_plot.md),
[`rp_size()`](rp_size.md)

## Author

Fred Hasselman
