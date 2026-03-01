# Inter-layer mutual information

Inter-layer mutual information

## Usage

``` r
mi_interlayer(g0, g1, probTable = FALSE)
```

## Arguments

- g0:

  An igraph object representing a layer in a multiplex graph

- g1:

  An igraph object representing a layer in a multiplex graph

- probTable:

  Option to return the table with marginal and joint degree distribution
  probabilities (default = `TRUE`)

## Value

The inter-layer mutual information between `g1` and `g2`. If
`probTable=TRUE`, a list object with two fields, the inter-layer mutual
information and the table with marginal and joint degree distributions

## Note

If the networks are weighted the strength distribution will be used
instead of the the degree distribution.

## See also

Other Redundancy measures (mutual information): [`mi_mat()`](mi_mat.md),
[`mif()`](mif.md)
