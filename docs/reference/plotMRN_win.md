# Plot windowed Multiplex Recurrence Network measures

Plot windowed Multiplex Recurrence Network measures

## Usage

``` r
plotMRN_win(
  df_mrn,
  layerMeasures = c("InterLayerMI", "InterLayerCor", "EdgeOverlap", "JRP")[1],
  vertexMeasures = c("none", "degree", "strength", "transitivity", "closeness",
    "betweenness")[1],
  weighted = FALSE,
  directed = FALSE,
  cumulative = FALSE,
  plotSD = FALSE,
  doPlot = TRUE
)
```

## Arguments

- df_mrn:

  Output from function [`mrn()`](mrn.md) with arguments set for a
  windowed analysis.

- layerMeasures:

  Character vector indicating which layer similarity measure(s) should
  be plotted. Valid elements in the vector are: `"InterLayerMI"`,
  `"InterLayerCor"`, `"EdgeOverlap"`, `"JRP"`.

- vertexMeasures:

  Character vector indicating which vertex measure(s) should be plotted.
  Valid elements in the vector are: `"none"`,
  `"degree"`,`"strength"`,`"transitivity"` (=
  clustering),`"closeness"`,`"betweenness"`. See the igraph::igraph The
  calculations depend on the value of `MRNweightedBy` that was used with
  function [`mrn()`](mrn.md).

- weighted:

  If `vertexMeasures` is not equal to `"none"`, should the graph be
  considered weighted? (default = `FALSE`)

- directed:

  If `vertexMeasures` is not equal to `"none"`, should the graph be
  considered directed? (default = `FALSE`)

- cumulative:

  If `vertexMeasures` is not equal to `"none"`, should the graph be
  considered cumulative? This will set the mode to `"upper"` (default =
  `FALSE`)

- plotSD:

  If `TRUE`, a ribbon with range mean ± SD will be plotted (default =
  `FALSE`)

- doPlot:

  Plot the igraph object. If `FALSE`, a just the graph object will be
  returned invisible.

## Value

a ggplot object

## See also

Other Tools for windowed analyses: [`ts_windower()`](ts_windower.md)
