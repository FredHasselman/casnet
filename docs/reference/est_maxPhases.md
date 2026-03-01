# Estimate the maximum number of Phases

Parameter sweep of function [rn_phases](rn_phases.md) for argument
`maxPhases`. Use to check at which value of `maxPhases` no additional
phases will be detected.

## Usage

``` r
est_maxPhases(
  RN,
  range = 2:10,
  minStatesinPhase = 1,
  maxStatesinPhase = NROW(RN),
  useDegree = FALSE,
  inverseWeight = TRUE,
  cleanUp = TRUE,
  removeSingularities = TRUE
)
```

## Arguments

- RN:

  Recurrence matrix

- range:

  Two element vector with minimum and maximum `c(min,max)` number of
  phases to check.

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

- useDegree:

  By default, node strength will be used to determine the phases. Set to
  `TRUE` to use the node degree (default = `FALSE`)

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

- cleanUp:

  Try to assign states to phases that were not assigned by the
  algorithm. If `FALSE`, these states will be added to category "Other"
  (default = `TRUE`)

- removeSingularities:

  Will remove states that recur only once (nodes with `degree(g) == 1`)
  (default = `FALSE`)

## Value

Data frame with maxPhases by detectedPhases
