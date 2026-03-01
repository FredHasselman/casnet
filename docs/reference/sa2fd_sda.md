# Informed Dimension estimate from SDA slope.

Conversion formula: Standardised Dispersion Analysis (SDA) estimate of
self-affinity parameter (`SA`) to an informed estimate of the fractal
dimension (FD).

## Usage

``` r
sa2fd_sda(sa, ...)
```

## Arguments

- sa:

  Self-afinity parameter estimate based on SDA slope (e.g.,
  [`fd_sda()`](fd_sda.md))).

- ...:

  Other arguments

## Value

An informed estimate of the Fractal Dimension, see Hasselman(2013) for
details.

## Details

Note that for some signals different PSD slope values project to a
single SDA slope. That is, SDA cannot distinguish dplyr::between all
variaties of power-law scaling in the frequency domain.

## References

Hasselman, F. (2013). When the blind curve is finite: dimension
estimation and model inference based on empirical waveforms. Frontiers
in Physiology, 4, 75. https://doi.org/10.3389/fphys.2013.00075

## Author

Fred Hasselman
