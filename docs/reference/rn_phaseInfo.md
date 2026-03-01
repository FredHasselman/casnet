# Extract Phases from weighted RN

This function will extract phases (regions of attraction in phase space)
based on a weighted `RN` object created with function [rn](rn.md). The
assumption is that coordinates in state space that are either close in
terms of distance, or are re-visited with short recurrence times, or
with high-frequency, are regions of attraction for the system.

## Usage

``` r
rn_phaseInfo(
  RN,
  maxPhases = NA,
  minStatesinPhase = 2,
  maxStatesinPhase = NROW(RN),
  selectionMethod = c("degree", "strength", "closeness", "betweenness")[1],
  inverseWeight = TRUE,
  returnCentroid = c("no", "mean.sd", "median.mad", "centroid")[1],
  removeSingularities = FALSE,
  standardise = c("none", "mean.sd", "median.mad", "unit")[4],
  returnGraph = FALSE,
  doPhaseProfilePlot = FALSE,
  plotCentroid = FALSE,
  dimNames = NULL,
  colOrder = FALSE,
  phaseColours = NULL,
  doSpiralPlot = FALSE,
  doPhaseSeriesPlot = FALSE,
  doPhaseDensityPlot = FALSE,
  doPhaseSpacePojectionPlot = FALSE,
  showPhaseSize = TRUE,
  showEpochLegend = TRUE,
  epochColours = NULL,
  epochLabel = "Phase",
  excludeVars = "",
  excludePhases = "",
  excludeTransients = FALSE,
  excludePhaseNeighbours = FALSE,
  excludeSingularities = FALSE,
  excludeNonrecurring = TRUE,
  alphaDensity = 0.4,
  splitFacets = NA,
  silent = FALSE
)
```

## Arguments

- RN:

  A matrix produced by the function [rn](rn.md)

- maxPhases:

  The maximum number of phases to extract. These will be the phases
  associated with the highest node degree or node strength. All other
  recurrent points will be labelled with "Other". If `NA`, the value
  will be set to `NROW(RN)`, this will return all the potential phases
  in the data irrespective of their frequency/strength of recurrence
  (default = `NROW(RN)`)

- minStatesinPhase:

  A parameter applied after the extraction of phases (limited by
  `maxPhases`). If any extracted phases do not have a minimum number of
  `minStatesinPhase` + `1` (= the state that was selected based on node
  strength), the phase will be removed from the result (default = `1`)

- maxStatesinPhase:

  A parameter applied after the extraction of phases (limited by
  `maxPhases`). If any extracted phases exceeds a maximum number of
  `maxStatesinPhase` + `1` (= the state that was selected based on node
  strength), the phase will be removed from the result (default =
  `NROW(RN)`)

- selectionMethod:

  How will the most "important" node be selected? Should be the name of
  an igraph::igraph function that returns local (vertex) measures
  (default = `"degree"`)

- inverseWeight:

  Whether to perform the operation `1/weight` on the edge weights. The
  default is `TRUE`, if the matrix was weighted by a distance metric
  (`weightedBy = "si"`) edges with smaller distances (recurring
  coordinates closer to the current coordinate) have greater impact on
  the node strength calculation used to select the phases. If the matrix
  was weighted by recurrence time (`weightedBy = "rt"`) and
  `inverseWeight = TRUE`, recurrent points with shorter recurrence times
  will have greater impact on the strength calculation. If
  `weightedBy = "rf"`, lower frequencies will end up having more impact.
  (default = `TRUE`)

- returnCentroid:

  Values can be `"no"`, `"mean.sd"`, `"median.mad"`, `"centroid"`. Any
  other value than `"no"` will return a data frame with the central
  tendency and deviation measures for each phase (default = `"no"`)

- removeSingularities:

  Will remove states that recur only once (nodes with `degree(g) == 1`)
  (default = `FALSE`)

- standardise:

  Standardise the series using [`ts_standardise()`](ts_standardise.md)
  with `adjustN = FALSE` (default = "mean.sd")

- returnGraph:

  Returns all the graph object objects of the plots that have been
  produced (default = `FALSE`)

- doPhaseProfilePlot:

  Produce a profile plot of the extracted phases (default = `TRUE`)

- plotCentroid:

  Plot the centroid requested in `returnCentroid`? (default = `FALSE`)

- dimNames:

  A vector of titles to use for the dimensions. If `NULL` the values
  will be read from the attribute of `RN`.

- colOrder:

  Should the order of the dimensions reflect the order of the columns in
  the dataset? If `FALSE` the order will be based on the values of the
  dimensions observed in the first extracted phase (default = `FALSE`)

- phaseColours:

  Colours for the different phases in the phase plot. If `epochColours`
  also has a value, `phaseColours` will be used instead (default =
  `NULL`)

- doSpiralPlot:

  Produce a plot of the recurrence network with the nodes coloured by
  phases (default = `FALSE`)

- doPhaseSeriesPlot:

  Produce a time series of the phases as they occur with a marginal
  histogram of their frequency (default = `FALSE`)

- showPhaseSize:

  Show the number of states in each phase in the labels (default = TRUE)

- showEpochLegend:

  Should a legend be shown for the epoch colours? (default = `TRUE`)

- epochColours:

  A vector of length `vcount(g)` with colour codes (default = `NULL`)

- epochLabel:

  A title for the epoch legend (default = `"Epoch"`)

- excludeVars:

  Exclude specific dimension variables by name. Leave empty to include
  all variables (default = `""`)

- excludePhases:

  Exclude Phases by their name (variable `phase_name`). Leave empty to
  include all Phases (after the other exclusion arguments) (default =
  `""`)

- excludeTransients:

  Should the category "Transient" be excluded from plots? (default =
  `FALSE`)

- excludePhaseNeighbours:

  Should the category "PhaseN" be excluded from plots? (default =
  `FALSE`)

- excludeSingularities:

  Should the category "Singularity" be excluded from plots? (default =
  `TRUE`)

- excludeNonrecurring:

  Should the category "Nonrecurring" be excluded from plots? (default =
  `TRUE`)

- alphaDensity:

  Alpha value for the density plots

- splitFacets:

  Integer value to indicate if sets of phases should be displayed in
  different facets? Pass an integer equalt to, or larger than the number
  of phases found to display all phases for which a density can be
  calculated (default = `NA`)

- silent:

  Silent-ish mode (default = `FALSE`)

- doPhaseSpaceProjectionPlot:

  produce a 2D `umap` projection of the phases (default = `FALSE`)

## Value

A data frame with information about the phases or a list object with
data and graph objects (if requested). The data frame contains the phase
name, number and size, as well as properties of the prototypical state
that was selected by the selection method, the node number/time
(`maxState_time`), degree (`maxState_degree`), strength
(`maxState_strength`).

## Details

The method used for the identification of phases is on the properties of
the `RN` object:

- If weighted by distance `"si"`, the inverse distance will be used,
  which means higher weights correspond to closer states.

- If weighted by recurrence time `"rt"`, the inverse time will be used,
  which means higher weights correspond to faster recurrence times.

The procedure is as follows:

1.  Identify the node with the highest strength

2.  Identify the nodes that connect to this node

3.  Identify the node with highest strength that does not connect to the
    node identified in step 1.

4.  Repeat until criteria set in `maxPhases`, `minStatesinPhase` and
    `maxStatesinPhase` are triggered.

## See also

Other Distance matrix operations (recurrence network):
[`mat_di2bi()`](mat_di2bi.md), [`mat_di2ch()`](mat_di2ch.md),
[`mat_di2we()`](mat_di2we.md), [`rn()`](rn.md),
[`rn_phases()`](rn_phases.md), [`rn_plot()`](rn_plot.md),
[`rn_recSpec()`](rn_recSpec.md)

## Examples

``` r
# Use the ManyAnalysts dataset to create a phase plot with default settings
data("manyAnalystsESM")
df <- manyAnalystsESM[4:10]
RN <- rn(y1 = df, doEmbed = FALSE, weighted = TRUE, weightedBy = "si", emRad = NA)
#> `emRad` was set to NA due to the value `weighted = TRUE`, if you want an unthresholded matrix set `weighted = FALSE`

# This returns 6 phases which have minimally 2 states
rn_phaseInfo(RN, doPhaseProfilePlot = TRUE)
#> 
#> ~~~o~~o~~casnet~~o~~o~~~
#> 
#> Recurring states with high similarity will be considered a phase
#> 
#> 
#> Looking for phases...
#> State at time 76 is template for phase 1 
#> State at time 1 is template for phase 2 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=9 [Phase 01.5] | State at t=37 [Phase 01.9] | State at t=40 [Phase 01.10] | State at t=48 [Phase 01.15] | State at t=122 [Phase 01.20] >> will be labelled as Transient(s)
#> 
#> State at time 17 is template for phase 3 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=4 [Phase 01.2] | State at t=10 [Phase 01.6] | State at t=13 [Phase 01.7] | State at t=35 [Phase 01.8] | State at t=46 [Phase 01.13] | State at t=47 [Phase 01.14] | State at t=122 [Phase 01.20] >> will be labelled as Transient(s)
#> 
#> State at time 75 is template for phase 4 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=100 [Phase 01.18] | State at t=119 [Phase 01.19] >> will be labelled as Transient(s)
#> 
#> State at time 38 is template for phase 5 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=2 [Phase 02.2] | State at t=36 [Phase 02.5] | State at t=55 [Phase 02.9] >> will be labelled as Transient(s)
#> 
#> State at time 54 is template for phase 6 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=2 [Phase 02.2] | State at t=43 [Phase 02.6] | State at t=55 [Phase 02.9] | State at t=107 [Phase 02.11] >> will be labelled as Transient(s)
#> 
#> State at time 85 is template for phase 7 
#> State at time 7 is template for phase 8 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=60 [Phase 07.3] >> will be labelled as Transient(s)
#> 
#> State at time 50 is template for phase 9 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=43 [Phase 02.6] | State at t=45 [Phase 02.7] | State at t=52 [Phase 02.8] | State at t=67 [Phase 02.10] >> will be labelled as Transient(s)
#> 
#> State at time 59 is template for phase 10 
#> State at time 39 is template for phase 11 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=4 [Phase 01.2] | State at t=13 [Phase 01.7] | State at t=47 [Phase 01.14] >> will be labelled as Transient(s)
#> 
#> State at time 34 is template for phase 12 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=20 [Phase 01.18] | State at t=36 [Phase 02.5] | State at t=100 [Phase 07.1] >> will be labelled as Transient(s)
#> 
#> State at time 89 is template for phase 13 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=16 [Phase 02.3] | State at t=43 [Phase 02.6] | State at t=67 [Phase 02.10] >> will be labelled as Transient(s)
#> 
#> State at time 106 is template for phase 14 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=79 [Phase 04.2] | State at t=105 [Phase 04.5] | State at t=116 [Phase 04.6] >> will be labelled as Transient(s)
#> 
#> State at time 51 is template for phase 15 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=2 [Phase 02.2] | State at t=117 [Phase 06.3] >> will be labelled as Transient(s)
#> 
#> State at time 78 is template for phase 16 
#> State at time 115 is template for phase 17 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=105 [Phase 04.5] >> will be labelled as Transient(s)
#> 
#> 
#> Found 17 phases with at least 2 states.


# Use min. number of states as the extraction criterion
rn_phaseInfo(RN, maxPhases = NA, minStatesinPhase = 7, doPhaseProfilePlot = TRUE)
#> 
#> ~~~o~~o~~casnet~~o~~o~~~
#> 
#> Recurring states with high similarity will be considered a phase
#> 
#> Warning: Phases with less than 7 states will be considered singularities.
#> 
#> Looking for phases...
#> State at time 76 is template for phase 1 
#> State at time 1 is template for phase 2 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=9 [Phase 01.5] | State at t=37 [Phase 01.9] | State at t=40 [Phase 01.10] | State at t=48 [Phase 01.15] | State at t=122 [Phase 01.20] >> will be labelled as Transient(s)
#> 
#> State at time 17 is template for phase 3 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=4 [Phase 01.2] | State at t=10 [Phase 01.6] | State at t=13 [Phase 01.7] | State at t=35 [Phase 01.8] | State at t=46 [Phase 01.13] | State at t=47 [Phase 01.14] | State at t=122 [Phase 01.20] >> will be labelled as Transient(s)
#> 
#> State at time 75 is template for phase 4 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=100 [Phase 01.18] | State at t=119 [Phase 01.19] >> will be labelled as Transient(s)
#> 
#> 
#> Found 4 phases with at least 7 states.

```
