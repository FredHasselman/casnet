# Estimate embedding lag (tau)

A wrapper for
[nonlinearTseries::timeLag](https://rdrr.io/pkg/nonlinearTseries/man/timeLag.html)

## Usage

``` r
est_emLag(y, selection.methods = "first.minimum", maxLag = length(y)/4, ...)
```

## Arguments

- y:

  Time series or numeric vector

- selection.methods:

  Selecting an optimal embedding lag (default: Return "first.e.decay",
  "first.zero", "first.minimum", "first.value", where value is 1/exp(1))

- maxLag:

  Maximal lag to consider (default: 1/4 of timeseries length)

- ...:

  Additional parameters

## Value

The ami function with requested minima

## See also

Other Estimate Recurrence Parameters: [`est_emDim()`](est_emDim.md),
[`est_parameters()`](est_parameters.md),
[`est_parameters_roc()`](est_parameters_roc.md),
[`est_radius()`](est_radius.md), [`est_radius_rqa()`](est_radius_rqa.md)
