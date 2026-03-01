# Estimate number of embedding dimensions

A wrapper for nonlinearTseries::estimateEmbeddingDim

## Usage

``` r
est_emDim(
  y,
  delay = est_emLag(y),
  maxDim = 15,
  threshold = 0.95,
  max.relative.change = 0.1,
  doPlot = FALSE,
  ...
)
```

## Arguments

- y:

  Time series or numeric vector

- delay:

  Embedding lag

- maxDim:

  Maximum number of embedding dimensions

- threshold:

  See
  [`nonlinearTseries::estimateEmbeddingDim()`](https://rdrr.io/pkg/nonlinearTseries/man/estimateEmbeddingDim.html)

- max.relative.change:

  See
  [`nonlinearTseries::estimateEmbeddingDim()`](https://rdrr.io/pkg/nonlinearTseries/man/estimateEmbeddingDim.html)

- doPlot:

  Plot

- ...:

  Other arguments (not in use)

## Value

Embedding dimensions

## Details

A wrapper for
[nonlinearTseries::estimateEmbeddingDim](https://rdrr.io/pkg/nonlinearTseries/man/estimateEmbeddingDim.html)

## See also

Other Estimate Recurrence Parameters: [`est_emLag()`](est_emLag.md),
[`est_parameters()`](est_parameters.md),
[`est_parameters_roc()`](est_parameters_roc.md),
[`est_radius()`](est_radius.md), [`est_radius_rqa()`](est_radius_rqa.md)
