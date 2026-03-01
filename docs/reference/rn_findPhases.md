# Find Phases

Algorithm to find similar recurring states that can be considered
phases, based on the node importance (`selectionMethod`) of a recurrence
network. This makes most sense if the state space dimension can be
interpreted.

## Usage

``` r
rn_findPhases(
  RN,
  weighted = NA,
  inverseWeight = TRUE,
  directed = NA,
  selectionMethod = c("degree", "strength", "closeness", "betweenness")[1],
  minStatesinPhase = 2,
  maxStatesinPhase = NROW(RN),
  maxPhases = NROW(RN),
  silent = FALSE
)
```

## Arguments

- RN:

  A matrix produced by the function [rn](rn.md)

- weighted:

  Should the matrix be considered to represent a weighted network?
  (default = `FALSE`)

- inverseWeight:

  Whether to perform the operation `1/weight` on the edge weights. This
  will only have an effect if wieght matters for the selection of most
  important nodes (e.g., `selectionMethod = "strength"`). The default is
  `TRUE`, so if the matrix was weighted by a distance metric
  (`weightedBy = "si"`) edges with smaller distances (recurring
  coordinates closer to the current coordinate) have greater impact on
  node strength. If the matrix was weighted by recurrence time
  (`weightedBy = "rt"`) and `inverseWeight = TRUE`, recurrent points
  with shorter recurrence times will have greater impact on the strength
  calculation. If `weightedBy = "rf"`, lower frequencies will end up
  having more impact. (default = `TRUE`)

- directed:

  Should the matrix be considered to represent a directed network?
  (default = `FALSE`)

- selectionMethod:

  How will the most "important" node be selected? Should be the name of
  an igraph::igraph function that returns local (vertex) measures
  (default = `"degree"`)

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

- maxPhases:

  The maximum number of phases to extract. These will be the phases
  associated with the highest node degree or node strength. All other
  recurrent points will be labelled with "Other". If `NA`, the value
  will be set to `NROW(RN)`, this will return all the potential phases
  in the data irrespective of their frequency/strength of recurrence
  (default = `10`)

- silent:

  Silent-ish mode

## Value

A data frame with nodes and the phase sequence.
