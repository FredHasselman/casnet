# Mutual Information variations

Mutual Information variations

## Usage

``` r
mi_mat(y, ID1, ID2, discreteBins = ceiling(2 * NROW(ID1)^(1/3)))
```

## Arguments

- y:

  A matrix with time series in columns

- ID1:

  ids

- ID2:

  ids

- discreteBins:

  Number of bins to use when discretizing the time series

## Value

mi in nats

## See also

Other Redundancy measures (mutual information):
[`mi_interlayer()`](mi_interlayer.md), [`mif()`](mif.md)
