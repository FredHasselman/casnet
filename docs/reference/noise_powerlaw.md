# Generate noise series with power law scaling exponent

Generate noise series with power law scaling exponent

## Usage

``` r
noise_powerlaw(
  y = NULL,
  alpha = -1,
  N = 512,
  standardise = FALSE,
  randomPower = FALSE,
  seed = NA
)
```

## Arguments

- y:

  Time series to use as a 'model'. If specified, `N` will be
  `N = length(y)`, and the series will be constructed based on
  `stats::fft(y)`.

- alpha:

  The log-log spectral slope, the scaling exponent. Use `0` for white
  noise, negative numbers for anti-persistant noises: `-1` for
  \\\frac{1}{f}\\ noise, positive numbers for persistent noises, e.g.
  `1` for blue noise.

- N:

  Length of the time series

- standardise:

  Forces scaling of the output to the range `[-1, 1]`, consequently the
  power law will not necessarily extend right down to `0Hz`.

- randomPower:

  If `TRUE` phases will be deterministic, uniformly distributed in
  `[-pi,pi]`. If `FALSE`, the spectrum will be stochastic with a
  Chi-square distribution. If `y` is not `NULL` this argument will be
  ignored.

- seed:

  Provide an integer number to set the seed for the random number
  generator in order to get reproducible results. If `NA` (default) no
  user defined seed will be set,

## Value

Time series with a power law of alpha.

## Note

This R code was adapted from a Matlab script called `powernoise.m` by
Max Little. The script contained the following commented text:

With no option strings specified, the power spectrum is
