# Plot Network Based on RQA

Plot Network Based on RQA

## Usage

``` r
plotNET_prep(
  g,
  labels = NA,
  nodeSize = c("degree", "hubscore", "strength", "eccentricity", "coreness")[1],
  rescaleSize = c(3, 10),
  nodeColour = TRUE,
  labelSize = "asnodesize",
  edgeWeight = "weight",
  edgeColour = FALSE,
  removeZeroDegree = TRUE,
  removeSelfLoops = TRUE,
  doPlot = TRUE
)
```

## Arguments

- g:

  An igraph object

- labels:

  Vertex labels. If `NA` is passed then first `V(g)$name` will be
  checked, for node labels. To create a plot without labels pass `NULL`
  (default = `NA`)

- nodeSize:

  Set node sizes by `degree(g, normalised = TRUE)` (default),
  `hubscore(g)$vector`, or, `strength(g)`, `eccentricity(g)`,
  `coreness(g)`. If a numeric value is passed all vertex sizes will be
  set to that value.

- rescaleSize:

  Use to rescale the measure indicated under `nodeSize` to `c(min,max)`
  for better node visibility. Pass `c(1,1)` for no rescaling (default =
  `c(3,10)`)

- nodeColour:

  Set to `TRUE` to colour nodes using the size of the node (default =
  `TRUE`)

- labelSize:

  Set labelsize: `"asnodesize"` sets the `cex` for the labels to
  coincide with nodesize (with min of .4 and max of 1.1). A single
  numeric value sets the `cex` of all labels to that value. A numeric
  vector of length two, `c(min,max)` will scale the label sizes to `min`
  and `max`

- edgeWeight:

  Set size of edges to `"E(g)$weight"` by passing "weight". If a single
  numeric value is provided all edges will be set to that value.

- edgeColour:

  Set to `TRUE` to colour edges using the weight of the edge (default =
  `FALSE`)

- removeZeroDegree:

  Remove vertices with `degree(g) == 0` (default = `TRUE`)

- removeSelfLoops:

  Calls `simplify(g)` (default = `TRUE`)

- doPlot:

  Plot the igraph object.

## Value

an igraph::igraph object

## See also

Other tools for plotting networks: [`plotNET_BA()`](plotNET_BA.md),
[`plotNET_SW()`](plotNET_SW.md),
[`plotNET_groupColour()`](plotNET_groupColour.md),
[`plotNET_groupWeight()`](plotNET_groupWeight.md)
