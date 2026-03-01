# Phase Density for each dimension

Create a ridge plot (requires package
[ggridges::ggridges](https://wilkelab.org/ggridges/reference/ggridges-package.html))

## Usage

``` r
plotRN_phaseDensities(
  phaseOutput,
  plotCentroid = NA,
  showEpochLegend = TRUE,
  epochColours = NULL,
  epochLabel = "Phase",
  excludeTransients = FALSE,
  excludePhaseNeighbours = FALSE,
  excludeSingularities = TRUE,
  excludeNonrecurring = TRUE,
  excludeVars = "",
  excludePhases = "",
  alphaDensity = 0.4,
  splitFacets = NA,
  showPhaseSize = TRUE,
  returnGraph = FALSE
)
```

## Arguments

- phaseOutput:

  Output from function [rn_phaseInfo](rn_phaseInfo.md)

- excludeTransients:

  Should the category "Transient" be excluded from plots? (default =
  `FALSE`)

- excludeSingularities:

  Should the category "Singularity" be excluded from plots? (default =
  `TRUE`)

- excludeNonrecurring:

  Should the category "Nonrecurring" be excluded from plots? (default =
  `TRUE`)

- excludeVars:

  Exclude specific dimension variables by name. Leave empty to include
  all variables (default = `""`)

- excludePhases:

  Exclude Phases by their name (variable `phase_name`). Leave empty to
  include all Phases (after the other exclusion arguments) (default =
  `""`)

- alphaDensity:

  Alpha value for the density plots

- splitFacets:

  Integer value to indicate if sets of phases should be displayed in
  different facets? Pass an integer equalt to, or larger than the number
  of phases found to display all phases for which a density can be
  calculated (default = `NA`)

- showPhaseSize:

  Show the number of states in each phase in the labels (default = TRUE)

- excludePhaseNeighbour:

  Should the category "PhaseNH" be excluded from plots? (default =
  `FALSE`)

## Value

data frame with phase series.

## Examples

``` r
RN <- rn(cumsum(rnorm(100)), emDim = 1, emLag = 1, emRad = NA, weighted = TRUE)
#> `emRad` was set to NA due to the value `weighted = TRUE`, if you want an unthresholded matrix set `weighted = FALSE`
#> Set `weightedBy` to 'si' due to the value `weighted = TRUE`
outPhases <- rn_phaseInfo(RN)
#> 
#> ~~~o~~o~~casnet~~o~~o~~~
#> 
#> Recurring states with high similarity will be considered a phase
#> 
#> 
#> Looking for phases...
#> State at time 4 is template for phase 1 
#> State at time 43 is template for phase 2 
#> State at time 18 is template for phase 3 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=10 [Phase 01.4] | State at t=11 [Phase 01.5] | State at t=19 [Phase 01.6] | State at t=21 [Phase 01.8] | State at t=50 [Phase 01.10] | State at t=66 [Phase 01.11] >> will be labelled as Transient(s)
#> 
#> State at time 49 is template for phase 4 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=2 [Phase 01.1] | State at t=5 [Phase 01.3] | State at t=20 [Phase 01.7] | State at t=39 [Phase 01.9] >> will be labelled as Transient(s)
#> 
#> State at time 100 is template for phase 5 
#> State at time 48 is template for phase 6 
#> State at time 91 is template for phase 7 
#> State at time 26 is template for phase 8 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=70 [Phase 05.1] >> will be labelled as Transient(s)
#> 
#> State at time 34 is template for phase 9 
#> State at time 28 is template for phase 10 
#> State at time 56 is template for phase 11 
#> State at time 80 is template for phase 12 
#> State at time 38 is template for phase 13 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=22 [Phase 04.1] | State at t=23 [Phase 04.2] | State at t=67 [Phase 04.5] >> will be labelled as Transient(s)
#> 
#> State at time 41 is template for phase 14 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=29 [Phase 10.2] | State at t=69 [Phase 10.4] >> will be labelled as Transient(s)
#> 
#> State at time 7 is template for phase 15 
#> State at time 27 is template for phase 16 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=73 [Phase 10.5] | State at t=87 [Phase 10.6] >> will be labelled as Transient(s)
#> 
#> State at time 35 is template for phase 17 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=40 [Phase 02.6] | State at t=47 [Phase 09.5] >> will be labelled as Transient(s)
#> 
#> 
#> Found 17 phases with at least 2 states.
plotRN_phaseDensities(outPhases)
#> Error in plotRN_phaseDensities(outPhases): object 'out' not found
```
