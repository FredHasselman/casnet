# Inter system recurrence networks

Inter system recurrence networks

## Usage

``` r
irn_plot(
  RMxy,
  gx = NULL,
  gy = NULL,
  merge_seed = 5,
  curvature = 0,
  angle = 0,
  scaleEdgeSize = 1,
  alphaE = 0.1,
  alphaV = 1,
  colorExy = "grey30",
  colorEyx = "grey70",
  sizeEbyCC = FALSE,
  doPlot = TRUE,
  ggplotReturn = FALSE,
  igraphReturn = FALSE
)
```

## Arguments

- RMxy:

  A cross recurrence matrix of x and y

- gx:

  An igraph object of recurrence network x

- gy:

  An igraph object of recurrence network y

- merge_seed:

  Seed

- curvature:

  The `curvature` parameter for edges see
  [`ggplot2::geom_curve()`](https://ggplot2.tidyverse.org/reference/geom_segment.html)
  (default = `-0.7`)

- angle:

  The `angle` parameter for edges see
  [`ggplot2::geom_curve()`](https://ggplot2.tidyverse.org/reference/geom_segment.html)
  (default = `90`)

- scaleEdgeSize:

  Scale the size of the edges by a constant:
  `E(g)$width * scaleEdgeSize` (default = `1/5`)

- alphaE:

  Set transparency for Edges. A single numeric, or a vector of length
  `ecount(g)` (default = `0.8`)

- alphaV:

  Set transparency for Vertices (default = `1`)

- colorExy:

  Edge colour x \>\> y

- colorEyx:

  Edge colour y \>\> x

- sizeEbyCC:

  Node size by Clustering Coefficient

- doPlot:

  Produce a plot? (default = `TRUE`)

- ggplotReturn:

  returns the ggplot object (default = `FALSE`)

- igraphReturn:

  returns the intermediate iGraph object. This will not look the same as
  the final graph, but has most of the attributes, like edge and vertex
  colors and spiral layout (default = `FALSE`)

## Value

A plot
