# Multiplex Recurrence Network

This function will create a Multiplex Recurrence Network from a list of
igraph::igraph objects that can be considered the layers of a network.
The layers must have the same number of nodes. There are two modes of
operation: *Layer similarity* (`MRNweightedBy` is set to
`"InterLayerMI"`, `"InterLayerCor"`, `"EdgeOvelap"`, or "`JRP`") and
*Layer importance* (not implemented yet). The former generates weighted
MRN based on *Interlayer Mutual Information*, *Interlayer Correlation*,
*Edge Overlap*, *Joint Recurrence Rate*.

## Usage

``` r
mrn(
  layers,
  MRNweightedBy = c("InterLayerMI", "InterLayerCor", "EdgeOverlap", "JRP")[1],
  jrp_measure = "RR",
  win = NA,
  step = NA,
  overlap = NA,
  alignment = "r",
  cumulative = FALSE,
  doPlot = FALSE,
  silent = TRUE
)
```

## Arguments

- layers:

  A list of igraph objects representing the layers of the multiplex
  network. The layer networks must all have the same number of vertices.

- MRNweightedBy:

  The measure to be used to evaluate the average structural similarities
  between the layers of the network. Valid options are: `"InterLayerMI"`
  (Mutual information based on similarity of the vertex degree across
  layers), `"InterLayerCor"` (pearson correlation between vertex degree
  or strength of each layer pair), `"EdgeOverlap"` (proportion of
  vertices sharing the same edges across layers), or `"JRP"` which will
  calculate a joint recurrence plot for each layer pair and add the
  value passed to `jrp_measure` as the measure of association (default =
  `InterLayerMI`)

- jrp_measure:

  Which JRP measure should be used to asses layer similarity. Use a
  column name of the data frame output by [rp_measures](rp_measures.md)
  (default = `"RR"`)

- win:

  The window size passed to [`ts_windower()`](ts_windower.md) in which
  to evaluate `"InterLayerMI"`, `"InterLayerCor"`, `"EdgeOvelap"`, or
  "`JRP`". (default = `NA`).

- step:

  The stepsize for the sliding window (default = `NA`).

- overlap:

  The window overlap passed to [`ts_windower()`](ts_windower.md) if
  `MRNweightedBy` is `"InterLayerMI"` or `"EdgeOvelap"`. The value of
  `step` will be ignored if `overlap` is not `NA`. (default = `NA`).

- alignment:

  Whether to right (`"r"`), center (`"c"`), or left (`"l"`) align the
  window.

- cumulative:

  To make the network represent cumulative time, set `directed = TRUE`
  and `cumulative = TRUE`. This will set the upper triangle of the
  recurrence matrix to `0` and ensures that the network edges represent
  recurrent values that have occurred in the `past` relative to the
  current observed value (node). If `directed = FALSE` the argument is
  ignored (default = `TRUE`).

- doPlot:

  Plot the multiplex recurrence network (default = `FALSE`).

- silent:

  Silent-ish mode. (default = FALSE).

## Value

A list object with a (windowed) Multiplex Recurrence Matrix based on
`MRNweightedBy`:

- *interlayerMI* - One or more matrices with edge weights between layers
  that represent the interlayer Mutual Information.

- *interlayerCor* - One or more matrices with edge weights between
  layers that represent the Pearson correlation between vertex degrees
  (or strength if layers are weighted graphs) of layers.

- *edgeOverlap* - One or more matrices with edge weights between layers
  that represent the overlapping edges between layers.

- *jointRP* - One or more matrices with edge weights between layers that
  represent the value of the joint RP passed in `jrp_measure`. The
  default is the `"RR"` of joint Recurrence Plot of the two layers.

- The field *meanValues* contains one or more matrices that represent
  the means and SDs of the weights requested in`MRNweightedBy`. The
  measure `eo_joint` (returned with `"edgeOverlap"`) refers to the
  number of edges shared among *all* layers of the MRN.

## Examples

``` r
# Create some layers
library(igraph)

layers <- list(g1 = igraph::sample_smallworld(1, 100, 5, 0.05),
g2 = igraph::sample_smallworld(1, 100, 5, 0.5),
g3 = igraph::sample_smallworld(1, 100, 5, 1))

mrn(layers = layers)

```
