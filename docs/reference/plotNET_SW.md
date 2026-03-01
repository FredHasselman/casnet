# Example of Strogatz-Watts small-world network

A wrapper around
[`igraph::sample_smallworld()`](https://r.igraph.org/reference/sample_smallworld.html)
with `dim=1`

## Usage

``` r
plotNET_SW(n = 100, k = 5, p = 0.05, doPlot = TRUE)
```

## Arguments

- n:

  Size of the lattice (integer)

- k:

  Neighbourhood size (integer)

- p:

  Rewiring probability (between `0` and `1`)

- doPlot:

  PLot the igraph object

## Value

A Strogatz-Watts small-world igraph object

## See also

[`igraph::sample_smallworld()`](https://r.igraph.org/reference/sample_smallworld.html)

Other tools for plotting networks: [`plotNET_BA()`](plotNET_BA.md),
[`plotNET_groupColour()`](plotNET_groupColour.md),
[`plotNET_groupWeight()`](plotNET_groupWeight.md),
[`plotNET_prep()`](plotNET_prep.md)
