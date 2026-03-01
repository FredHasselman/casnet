# Estimate Radius without building a recurrence matrix

Find a fixed radius without building the recurrence matrix.

## Usage

``` r
est_radius_rqa(
  y1 = NULL,
  y2 = NULL,
  emDim = NA,
  emLag = NA,
  AUTO = NULL,
  method = "Euclidean",
  startRadius = NULL,
  targetValue = 0.05,
  tol = 0.01,
  maxIter = 100,
  theiler = NA,
  histIter = FALSE,
  standardise = c("mean.sd", "median.mad", "none")[3],
  radiusOnFail = c("tiny", "huge", "percentile")[3],
  silent = FALSE,
  useParallel = TRUE,
  doEmbed = TRUE
)
```

## Arguments

- y1:

  A numeric vector or time series

- y2:

  A numeric vector or time series for cross recurrence

- emDim:

  The embedding dimensions

- emLag:

  The embedding lag

- AUTO:

  Auto-recurrence? (default = `FALSE`)

- method:

  Distance measure to use. Any option that is valid for argument
  `method` of `proxy::dist()`. Type `proxy::pr_DB$get_entries()` to see
  a list of all the options. Common methods are:
  `"Euclidean", "Manhattan", "Minkowski", "Chebysev"` (or the same but
  shorter: `"L2","L1","Lp", "max"` distance). To use the shape based
  distance for phase-based recurrence use `"SBD"` (default =
  `"Euclidean"`)

- startRadius:

  The starting value for the radius (default = SD of time series values)

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

- standardise:

  Standardise `y` if `type == "optimal"`

- radiusOnFail:

  Radius to return when search fails `"tiny" = 0 + ,Machine.double.eps`,
  this will likely cause a matrix full of zeros.
  `"huge" = 1 + max. distance`, which will give a matrix full of ones,
  `"minimum" = minimum distance in matrix`.

- silent:

  Silent-ish

- useParallel:

  Should evaluation run using package parallel? This is will only be
  beneficial if the time series contains more than 10k data points
  (default = `TRUE`)

- doEmbed:

  If `FALSE`, a distance matrix will be returned that is not embedded by
  `emDim` and `emLag` (Multidimensional RQA). If `y1` and/or `y2` are
  data frames, the columns will be used as the state space dimensions
  (default = `TRUE`)

- noiseLevel:

  Noise level to construct the `signal + noiseLevel *`
  \\N(\mu=0,\sigma=1)\\ (default = `0.75`)

- noiseType:

  Type

- plotROC:

  Generates an ROC plot if `type = "optimal"`

## Value

A dataframe listing settings used to search for the radius, the radius
found given the settings and the recurrence rate produced by the radius
(either 1 row or the entire iteration history)

## See also

Other Estimate Recurrence Parameters: [`est_emDim()`](est_emDim.md),
[`est_emLag()`](est_emLag.md), [`est_parameters()`](est_parameters.md),
[`est_parameters_roc()`](est_parameters_roc.md),
[`est_radius()`](est_radius.md)
