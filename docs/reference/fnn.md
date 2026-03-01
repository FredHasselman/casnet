# False Nearest Neighbours

Search for FNN to get an optimal Embedding Dimension using by using
[`nonlinearTseries::findAllNeighbours()`](https://rdrr.io/pkg/nonlinearTseries/man/findAllNeighbours.html)
in a loop.

## Usage

``` r
fnn(y, emLag = 1, maxDim = 10, radius = sd(y)/10, number.boxes = NULL)
```

## Arguments

- y:

  A numeric vector or time series

- emLag:

  Optimal embedding lag (delay), e.g., provided by an optimising
  algorithm. If `NULL` the lags based on the mutual information in
  `lagMethods` will be reported. If a numeric value representing a valid
  lag is passed, this value will be used to estimate the number of
  dimensions (default = `NULL`)

- maxDim:

  Maximum number of embedding dimensions to consider (default = `10`)

- radius:

  Size of the neighbourhood: Every point smaller than the radius will be
  considered a near neighbour, see
  [`nonlinearTseries::findAllNeighbours()`](https://rdrr.io/pkg/nonlinearTseries/man/findAllNeighbours.html)
  (default = `sd(y)/10`).

- number.boxes:

  Integer representing number of boxes to to speed up neighbour search,
  if `NULL` an optimal number will be chosen
  [`nonlinearTseries::findAllNeighbours()`](https://rdrr.io/pkg/nonlinearTseries/man/findAllNeighbours.html)
  (default = `NULL`).

## Value

FNN curve
