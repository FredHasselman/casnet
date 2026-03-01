# Relative Roughness

Relative Rougness is a ratio of local variance (autocovariance at lag-1)
to global variance (autocovariance at lag-0) that can be used to
classify different 'noises'.

## Usage

``` r
fd_RR(y)
```

## Arguments

- y:

  A numeric vector.

## Value

The Relative Roughness of y, the values of local and global variance are
returned as attributes

## Details

\$\$RR = 2 \* \left\[1 - \frac{\gamma(y)}{Var(y)}\right\]\$\$

## References

- Marmelat, V., Torre, K., & Delignieres, D. (2012). Relative roughness:
  an index for testing the suitability of the monofractal model.
  *Frontiers in Physiology, 3*, 208.

## See also

Other Fluctuation Analyses: [`fd_allan()`](fd_allan.md),
[`fd_dfa()`](fd_dfa.md), [`fd_mfdfa()`](fd_mfdfa.md),
[`fd_psd()`](fd_psd.md), [`fd_sda()`](fd_sda.md),
[`fd_sev()`](fd_sev.md)
