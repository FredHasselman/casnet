# Estimate RQA parameters

Find optimal parameters for constructing a Recurrence Matrix. A wrapper
for various algorithms used to find optimal values for the embedding
delay and the number of embedding dimensions.

## Usage

``` r
est_parameters(
  y,
  lagMethods = c("first.minimum", "global.minimum", "max.lag"),
  estimateDimensions = "preferSmallestInLargestHood",
  maxDim = 10,
  emLag = NULL,
  maxLag = NA,
  minVecLength = 20,
  nnSize = NA,
  nnRadius = NA,
  nnThres = 10,
  theiler = 0,
  doPlot = TRUE,
  silent = TRUE,
  ...
)
```

## Arguments

- y:

  A numeric vector or time series

- lagMethods:

  A character vector with one or more of the following strings:
  `"first.minimum","global.minimum","max.lag"`. If `emLag` represents a
  valid lag this value will be reported as `"user.lag"` (default =
  `c("first.minimum","global.minimum","max.lag")`)

- estimateDimensions:

  Decide on an optimal embedding dimension relative to the values in
  `maxDim` and `lagMethods`, according to a number of preferences passed
  as a character vector. The order in which the preferences appear in
  the vector affects the selection procedure, with index `1` being most
  important preference. The following options are available:

  - `preferNone` - No optimal number will be picked all other
    preferences will be ignored

  - `preferSmallestDim` - Pick smallest number of dimensions associated
    with a percentage NN below `nnThres`

  - `preferSmallestNN` - Pick the number of dimensions that is
    associated with the smallest percentage NN below `nnThres`

  - `preferSmallestLag` - If the value of `nnThres` does not lead to a
    unique preference for a pair of dimension and lag values, use the
    pair with the smallest lag

  - `preferSmallestInLargestHood` - The default option: If no unique
    pair can be found, prefer pairs with smallest values for lag,
    dimensions, percentage NN for the largest NN size

- maxDim:

  Maximum number of embedding dimensions to consider (default = `10`)

- emLag:

  Optimal embedding lag (delay), e.g., provided by an optimising
  algorithm. If `NULL` the lags based on the mutual information in
  `lagMethods` will be reported. If a numeric value representing a valid
  lag is passed, this value will be used to estimate the number of
  dimensions (default = `NULL`)

- maxLag:

  Maximum embedding lag to consider. If `NA` then the value is
  caclulated as `floor(length(y)/(maxDim+1))` (default = `NA`)

- minVecLength:

  The minimum length of state space vectors after delay-embedding. For
  short time series, this will affect the possible values of `maxDim`
  that can be used to evaluate the drop in nearest neighbours. In
  general it is not recommended to evaluate high dimensional state
  spaces, based on a small number of state soace coordinates, the
  default is an absolute minimum and possibly even lower than that.
  (default = `20`)

- nnSize:

  Neighbourhood diameter (integer, the `number.boxes` parameter of
  [`tseriesChaos::false.nearest()`](https://rdrr.io/pkg/tseriesChaos/man/false.nearest.html))
  used to speed up neighbour search. (default = `NA`)

- nnRadius:

  Points smaller than the radius are considered neighbours. If `NA` the
  value will be `sd(y)/10` (default = `NA`)

- nnThres:

  Threshold value (in percentage 0-100) representing the percentage of
  Nearest Neighbours that would be acceptable when using N surrogate
  dimensions. The smallest number of surrogate dimensions that yield a
  value below the threshold will be considered optimal (default = `10`)

- theiler:

  Theiler window on distance matrix (default = `0`)

- doPlot:

  Produce a diagnostic plot the results (default = `TRUE`)

- silent:

  Silent-ish mode

- ...:

  Other parameters passed to
  [`nonlinearTseries::timeLag()`](https://rdrr.io/pkg/nonlinearTseries/man/timeLag.html)

## Value

A list object containing the optimal values (as indicated by the user)
and iteration history.

## Details

A number of functions are called to determine optimal parameters for
delay embedding a time series:

- Embedding lag (`emLag`): The default is to call
  [`est_emLag()`](est_emLag.md), which is a wrapper around
  [`nonlinearTseries::timeLag()`](https://rdrr.io/pkg/nonlinearTseries/man/timeLag.html)
  with `technique=ami` to get lags based on the mutual information
  function.

- Embedding dimension (`m`, `emDim`): The default is to call
  [`est_emDim()`](est_emDim.md), which is a wrapper around
  [`tseriesChaos::false.nearest()`](https://rdrr.io/pkg/tseriesChaos/man/false.nearest.html)

## See also

Other Estimate Recurrence Parameters: [`est_emDim()`](est_emDim.md),
[`est_emLag()`](est_emLag.md),
[`est_parameters_roc()`](est_parameters_roc.md),
[`est_radius()`](est_radius.md), [`est_radius_rqa()`](est_radius_rqa.md)

## Examples

``` r
set.seed(4321)
est_parameters(y=rnorm(100))
#> Error in which(NROW(y) - (1:maxDim * maxLag) %>=% minVecLength): argument to 'which' is not logical
```
