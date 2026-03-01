# Example of Barabasi scale-free network

A wrapper around
[`igraph::sample_pa()`](https://r.igraph.org/reference/sample_pa.html)

## Usage

``` r
plotNET_BA(n = 100, pwr = 1, out.dist = NULL, doPlot = TRUE)
```

## Arguments

- n:

  Number of vertices

- pwr:

  Power of preferential attachment

- out.dist:

  Degree distribution

- doPlot:

  Plot the igraph object

## Value

A Barabasi scale-free igraph object

## See also

[`igraph::sample_pa()`](https://r.igraph.org/reference/sample_pa.html)

Other tools for plotting networks: [`plotNET_SW()`](plotNET_SW.md),
[`plotNET_groupColour()`](plotNET_groupColour.md),
[`plotNET_groupWeight()`](plotNET_groupWeight.md),
[`plotNET_prep()`](plotNET_prep.md)
