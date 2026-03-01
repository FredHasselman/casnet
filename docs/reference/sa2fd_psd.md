# Informed Dimension estimate from Spectral Slope (aplha)

Conversion formula: From periodogram based self-affinity parameter
estimate (`sa`) to an informed estimate of the (fractal) dimension (FD).

## Usage

``` r
sa2fd_psd(sa, ...)
```

## Arguments

- sa:

  Self-Affinity parameter estimate based on PSD slope (e.g.,
  [`fd_psd()`](fd_psd.md)))

- ...:

  Other arguments

## Value

An informed estimate of the Fractal Dimension, see Hasselman(2013) for
details.

## Details

The spectral slope will be converted to a dimension estimate using:

\$\$D\_{PSD}\approx\frac{3}{2}+\frac{14}{33}\*\tanh\left(Slope \*
\ln(1+\sqrt{2})\right)\$\$

## References

Hasselman, F. (2013). When the blind curve is finite: dimension
estimation and model inference based on empirical waveforms. Frontiers
in Physiology, 4, 75. https://doi.org/10.3389/fphys.2013.00075

## Author

Fred Hasselman
