# Estimate Radius.

Find a fixed or optimal radius.

## Usage

``` r
est_radius(
  RM = NULL,
  y1 = NULL,
  y2 = NULL,
  emLag = NA,
  emDim = NA,
  doEmbed = TRUE,
  method = "Euclidean",
  type = c("fixed", "optimal")[1],
  startRadius = NULL,
  eachRadius = 1,
  targetMeasure = c("RR", "DET", "LAM", "T1", "all")[1],
  targetValue = 0.05,
  tol = 0.01,
  maxIter = 100,
  theiler = NA,
  histIter = FALSE,
  noiseLevel = 0.75,
  noiseType = c("normal", "uniform")[1],
  plotROC = FALSE,
  standardise = c("mean.sd", "median.mad", "none")[3],
  radiusOnFail = c("tiny", "huge", "percentile")[1],
  silent = FALSE
)
```

## Arguments

- RM:

  Unthresholded Recurrence Matrix

- y1:

  A numeric vector or time series

- y2:

  A numeric vector or time series for cross recurrence

- emLag:

  The embedding lag

- emDim:

  The embedding dimensions

- doEmbed:

  If `FALSE`, a distance matrix will be returned that is not embedded by
  `emDim` and `emLag` (Multidimensional RQA). If `y1` and/or `y2` are
  data frames, the columns will be used as the state space dimensions
  (default = `TRUE`)

- method:

  Distance measure to use. Any option that is valid for argument
  `method` of `proxy::dist()`. Type `proxy::pr_DB$get_entries()` to see
  a list of all the options. Common methods are:
  `"Euclidean", "Manhattan", "Minkowski", "Chebysev"` (or the same but
  shorter: `"L2","L1","Lp", "max"` distance). To use the shape based
  distance for phase-based recurrence use `"SBD"` (default =
  `"Euclidean"`)

- type:

  Either `"fixed"` (default) or `"optimal"`, `"fixed"` will search for a
  radius that is close to the value for the `targetMeasure` in
  `targetValue`, `"optimal"` will optimise the radius for the
  `targetMeasure`, `targetValue` is ignored.

- startRadius:

  If `type = "fixed"` this is the starting value for the radius (default
  = percentile of unique distances in RM given by `targetValue`). If
  `type = "optimal"` this will be a range of radius values (in
  normalised SD units) that will be considered (default =
  `seq(0,2,by=.01)`)

- eachRadius:

  If `type = "optimal"` this is the number of signal and noise series
  that will be generated for each level in `startRadius` (default = `1`)

- targetMeasure:

  If `type = "optimal"`, it must be a character vector indicating which
  recurrence measure to optimise the radius for, options are "RR"
  (default), "DET", "LAM", "T1", and "all". The option
  `targetMeasure = "all"` will report all the optimal values obtained
  from one realisation of `startRadius * eachRadius` signal and noise
  series.

- targetValue:

  When argument `type` is set to "fixed", the value represents the
  target value for the measure in `targetMeasure` (default =
  `RR = .05`).

- tol:

  Tolerance for achieving `targetValue` for `targetMeasure` (default =
  `0.01`)

- maxIter:

  If `type = "fixed"`: Maximum number of iterations to reach
  targetValue.

- theiler:

  Size of theiler window (default `0`)

- histIter:

  Return iteration history? (default = `FALSE`)

- noiseLevel:

  Noise level to construct the `signal + noiseLevel *`
  \\N(\mu=0,\sigma=1)\\ (default = `0.75`)

- noiseType:

  Type

- plotROC:

  Generates an ROC plot if `type = "optimal"`

- standardise:

  Standardise `y` if `type == "optimal"`

- radiusOnFail:

  Radius to return when search fails `"tiny" = 0 + ,Machine.double.eps`,
  this will likely cause a matrix full of zeros.
  `"huge" = 1 + max. distance in RM`, which will give a matrix full of
  ones,
  `"percentile" = quantile(RM, prob = targetValue) of distances greater than 0`.

- silent:

  Silent-ish

## Value

A dataframe listing settings used to search for the radius, the radius
found given the settings and the recurrence rate produced by the radius
(either 1 row or the entire iteration history)

## See also

Other Estimate Recurrence Parameters: [`est_emDim()`](est_emDim.md),
[`est_emLag()`](est_emLag.md), [`est_parameters()`](est_parameters.md),
[`est_parameters_roc()`](est_parameters_roc.md),
[`est_radius_rqa()`](est_radius_rqa.md)
