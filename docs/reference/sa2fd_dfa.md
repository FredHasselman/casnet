# Informed Dimension estimate from DFA slope (H)

Conversion formula: Detrended Fluctuation Analysis (DFA) estimate of the
Hurst exponent (a self-affinity parameter `sa`) to an informed estimate
of the (fractal) dimension (FD).

## Usage

``` r
sa2fd_dfa(sa, ...)
```

## Arguments

- sa:

  Self-Afinity parameter estimate based on DFA slope (e.g.,
  [`fd_sda()`](fd_sda.md))).

- ...:

  Other arguments

## Value

An informed estimate of the Fractal Dimension, see Hasselman(2013) for
details.

## Details

The DFA slope (H) will be converted to a dimension estimate using:

\$\$D\_{DFA}\approx 2-(\tanh(\log(3)\*sa)) \$\$

## References

Hasselman, F. (2013). When the blind curve is finite: dimension
estimation and model inference based on empirical waveforms. Frontiers
in Physiology, 4, 75. https://doi.org/10.3389/fphys.2013.00075

## Author

Fred Hasselman
