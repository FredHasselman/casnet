# Vertex and Edge Group Colours

Identify Vertex and/or Edge groups by colour.

## Usage

``` r
plotNET_groupColour(
  g,
  groups = NULL,
  colourV = TRUE,
  alphaV = 1,
  colourE = FALSE,
  alphaE = 0.8,
  groupColours = NULL,
  defaultEdgeColour = "grey70",
  doPlot = TRUE,
  returnPlot = FALSE
)
```

## Arguments

- g:

  An igraph object

- groups:

  A named numeric vector with `length(V(g))` integers representing each
  group, or, a named character vector describing each group. If
  `names(groups)==NULL` then the names of the vector will be set as
  `names(groups) == V(g)$name`. If `V(g)$name==NULL`, the names of the
  vector will be set by the Vertex index

- colourV:

  Colour Vertices based on `groups` (default = `TRUE`)

- alphaV:

  Set transparency for Vertices (default = `1`)

- colourE:

  Colour Edges based on `groups`. Edges connecting to vertices of the
  same group will be coloured as the group (default = `FALSE`)

- alphaE:

  Set transparency for Edges. A single numeric, or a vector of length
  `ecount(g)` (default = `0.8`)

- groupColours:

  A list of length `groups` with valid colour codes

- defaultEdgeColour:

  Default edge colour

- doPlot:

  Plot the igraph object

- returnPlot:

  return the
  [ggplot2::ggplot2](https://ggplot2.tidyverse.org/reference/ggplot2-package.html)
  object

## Value

An igraph object with vertices and/or edges coloured by groups listed in
`groups`

## See also

Other tools for plotting networks: [`plotNET_BA()`](plotNET_BA.md),
[`plotNET_SW()`](plotNET_SW.md),
[`plotNET_groupWeight()`](plotNET_groupWeight.md),
[`plotNET_prep()`](plotNET_prep.md)
