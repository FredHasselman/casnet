# Profile Plot

Plot a profile (values of the dimensions) for each phase.

## Usage

``` r
plotRN_phaseProfiles(
  phaseOutput,
  plotCentroid = NA,
  showEpochLegend = TRUE,
  epochColours = NULL,
  epochLabel = "Phase",
  excludeTransients = FALSE,
  excludePhaseNeighbours = FALSE,
  excludeNonrecurring = TRUE,
  excludeSingularities = TRUE,
  excludeVars = "",
  excludePhases = "",
  showPhaseSize = TRUE,
  returnGraph = FALSE,
  colOrder = NA
)
```

## Arguments

- phaseOutput:

  Output from function [rn_phaseInfo](rn_phaseInfo.md)

- excludeTransients:

  Should the category "Transient" be excluded from plots? (default =
  `FALSE`)

- excludeNonrecurring:

  Should the category "Nonrecurring" be excluded from plots? (default =
  `TRUE`)

- excludeSingularities:

  Should the category "Singularity" be excluded from plots? (default =
  `TRUE`)

- excludeVars:

  Exclude specific dimension variables by name. Leave empty to include
  all variables (default = `""`)

- excludePhases:

  Exclude Phases by their name (variable `phase_name`). Leave empty to
  include all Phases (after the other exclusion arguments) (default =
  `""`)

- showPhaseSize:

  Show the number of states in each phase in the labels (default = TRUE)

## Value

plot

## Examples

``` r
RN <- rn(y1 = data.frame(x = rnorm(100), y= rnorm(100)), weighted = TRUE)
#> `emRad` was set to NA due to the value `weighted = TRUE`, if you want an unthresholded matrix set `weighted = FALSE`
#> Set `weightedBy` to 'si' due to the value `weighted = TRUE`
phase_out <- rn_phaseInfo(RN, returnCentroid = "mean.sd")
#> 
#> ~~~o~~o~~casnet~~o~~o~~~
#> 
#> Recurring states with high similarity will be considered a phase
#> 
#> 
#> Looking for phases...
#> State at time 28 is template for phase 1 
#> State at time 15 is template for phase 2 
#> State at time 41 is template for phase 3 
#> State at time 17 is template for phase 4 
#> State at time 91 is template for phase 5 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=9 [Phase 02.2] | State at t=82 [Phase 02.10] >> will be labelled as Transient(s)
#> 
#> State at time 65 is template for phase 6 
#> State at time 5 is template for phase 7 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=33 [Phase 01.5] | State at t=46 [Phase 01.9] | State at t=66 [Phase 01.12] >> will be labelled as Transient(s)
#> 
#> State at time 76 is template for phase 8 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=1 [Phase 01.1] | State at t=9 [Phase 01.14] | State at t=67 [Phase 02.2] | State at t=98 [Phase 05.4] >> will be labelled as Transient(s)
#> 
#> State at time 90 is template for phase 9 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=8 [Phase 03.1] | State at t=22 [Phase 03.2] | State at t=63 [Phase 03.10] >> will be labelled as Transient(s)
#> 
#> State at time 100 is template for phase 10 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=27 [Phase 06.2] | State at t=50 [Phase 06.3] | State at t=79 [Phase 06.6] >> will be labelled as Transient(s)
#> 
#> State at time 42 is template for phase 11 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=58 [Phase 03.7] | State at t=59 [Phase 03.8] | State at t=85 [Phase 03.11] >> will be labelled as Transient(s)
#> 
#> State at time 23 is template for phase 12 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=18 [Phase 04.3] | State at t=68 [Phase 04.6] >> will be labelled as Transient(s)
#> 
#> State at time 70 is template for phase 13 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=61 [Phase 10.1] >> will be labelled as Transient(s)
#> 
#> State at time 80 is template for phase 14 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=31 [Phase 04.5] | State at t=95 [Phase 04.8] >> will be labelled as Transient(s)
#> 
#> State at time 4 is template for phase 15 
#> State at time 69 is template for phase 16 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=66 [Phase 01.12] >> will be labelled as Transient(s)
#> 
#> State at time 89 is template for phase 17 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=73 [Phase 06.5] | State at t=87 [Phase 06.7] >> will be labelled as Transient(s)
#> 
#> 
#> Found 17 phases with at least 2 states.
plotRN_phaseProfiles(phaseOutput = phase_out)
#> Warning: Removed 2 rows containing missing values or values outside the scale range
#> (`geom_segment()`).
#> Warning: Removed 2 rows containing missing values or values outside the scale range
#> (`geom_segment()`).


```
