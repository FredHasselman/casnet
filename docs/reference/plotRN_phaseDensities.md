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
#> State at time 25 is template for phase 1 
#> State at time 27 is template for phase 2 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=48 [Phase 01.4] | State at t=49 [Phase 01.5] | State at t=51 [Phase 01.7] | State at t=68 [Phase 01.8] | State at t=79 [Phase 01.11] >> will be labelled as Transient(s)
#> 
#> State at time 52 is template for phase 3 
#> State at time 59 is template for phase 4 
#> State at time 14 is template for phase 5 
#> State at time 57 is template for phase 6 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=62 [Phase 03.3] | State at t=82 [Phase 03.9] >> will be labelled as Transient(s)
#> 
#> State at time 78 is template for phase 7 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=64 [Phase 03.4] | State at t=66 [Phase 03.5] | State at t=81 [Phase 03.8] >> will be labelled as Transient(s)
#> 
#> State at time 96 is template for phase 8 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=56 [Phase 04.4] | State at t=63 [Phase 04.5] | State at t=84 [Phase 04.7] | State at t=89 [Phase 06.2] | State at t=95 [Phase 06.5] >> will be labelled as Transient(s)
#> 
#> State at time 67 is template for phase 9 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=23 [Phase 01.2] | State at t=24 [Phase 01.6] | State at t=50 [Phase 01.9] | State at t=76 [Phase 07.2] | State at t=87 [Phase 07.4] >> will be labelled as Transient(s)
#> 
#> State at time 15 is template for phase 10 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=9 [Phase 05.1] | State at t=10 [Phase 05.2] | State at t=19 [Phase 05.5] | State at t=20 [Phase 05.6] >> will be labelled as Transient(s)
#> 
#> State at time 26 is template for phase 11 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=44 [Phase 05.7] | State at t=45 [Phase 05.8] >> will be labelled as Transient(s)
#> 
#> State at time 90 is template for phase 12 
#> State at time 12 is template for phase 13 
#> State at time 33 is template for phase 14 
#> State at time 53 is template for phase 15 
#> State at time 61 is template for phase 16 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=4 [Phase 04.1] | State at t=93 [Phase 15.4] >> will be labelled as Transient(s)
#> 
#> State at time 91 is template for phase 17 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=60 [Phase 12.5] | State at t=86 [Phase 15.2] | State at t=100 [Phase 15.3] >> will be labelled as Transient(s)
#> 
#> State at time 18 is template for phase 18 
#> State at time 22 is template for phase 19 
#> 
#> Found 19 phases with at least 2 states.
plotRN_phaseDensities(outPhases)
#> Error in plotRN_phaseDensities(outPhases): object 'out' not found
```
