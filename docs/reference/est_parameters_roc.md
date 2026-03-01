# Estimate ROC radius

Experimental.

## Usage

``` r
est_parameters_roc(
  y,
  emRad,
  emDim = 1,
  emLag = 1,
  noiseLevel = 0.75,
  standardise = c("mean.sd", "median.mad", "none")[3],
  noiseType = c("normal", "uniform")[1]
)
```

## Arguments

- y:

  y

- emRad:

  radius

- emDim:

  embedding Dims

- emLag:

  embedding Lag

- noiseLevel:

  noise Level

- standardise:

  Standardise y? Choose from "mean.sd","median.mad","none".

- noiseType:

  Use a Normal distribution of uniform distribution for noiselevels

## Value

data frame for ROC

## See also

Other Estimate Recurrence Parameters: [`est_emDim()`](est_emDim.md),
[`est_emLag()`](est_emLag.md), [`est_parameters()`](est_parameters.md),
[`est_radius()`](est_radius.md), [`est_radius_rqa()`](est_radius_rqa.md)
