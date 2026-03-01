# Phase Series plot

Plot the sequence of phases as a time series.

## Usage

``` r
plotRN_phaseTimeSeries(
  phaseOutput,
  showEpochLegend = TRUE,
  epochColours = NULL,
  epochLabel = "Phase",
  excludeVars = "",
  excludePhases = "",
  returnGraph = FALSE
)
```

## Arguments

- phaseOutput:

  Output from function [rn_phaseInfo](rn_phaseInfo.md)

- excludeVars:

  Exclude specific dimension variables by name. Leave empty to include
  all variables (default = `""`)

- excludePhases:

  Exclude Phases by their name (variable `phase_name`). Leave empty to
  include all Phases (after the other exclusion arguments) (default =
  `""`)

## Value

phase density plot

## Examples

``` r
RN <- rn(y1 = rnorm(100), weighted = TRUE)
#> `emRad` was set to NA due to the value `weighted = TRUE`, if you want an unthresholded matrix set `weighted = FALSE`
#> Set `weightedBy` to 'si' due to the value `weighted = TRUE`
phase_out <- rn_phaseInfo(RN)
#> 
#> ~~~o~~o~~casnet~~o~~o~~~
#> 
#> Recurring states with high similarity will be considered a phase
#> 
#> 
#> Looking for phases...
#> State at time 39 is template for phase 1 
#> State at time 17 is template for phase 2 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=22 [Phase 01.1] | State at t=43 [Phase 01.5] | State at t=63 [Phase 01.10] >> will be labelled as Transient(s)
#> 
#> State at time 12 is template for phase 3 
#> State at time 33 is template for phase 4 
#> State at time 42 is template for phase 5 
#> State at time 2 is template for phase 6 
#> State at time 28 is template for phase 7 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=37 [Phase 01.3] | State at t=44 [Phase 01.6] | State at t=52 [Phase 01.7] | State at t=59 [Phase 01.9] | State at t=97 [Phase 01.13] >> will be labelled as Transient(s)
#> 
#> State at time 86 is template for phase 8 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=3 [Phase 02.1] | State at t=48 [Phase 02.5] | State at t=67 [Phase 02.6] | State at t=84 [Phase 02.7] >> will be labelled as Transient(s)
#> 
#> State at time 1 is template for phase 9 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=15 [Phase 06.3] | State at t=40 [Phase 06.4] >> will be labelled as Transient(s)
#> 
#> State at time 70 is template for phase 10 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=14 [Phase 03.6] | State at t=26 [Phase 04.1] | State at t=68 [Phase 04.2] | State at t=73 [Phase 04.6] >> will be labelled as Transient(s)
#> 
#> State at time 4 is template for phase 11 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=20 [Phase 08.1] | State at t=100 [Phase 08.3] >> will be labelled as Transient(s)
#> 
#> State at time 38 is template for phase 12 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=98 [Phase 06.7] >> will be labelled as Transient(s)
#> 
#> State at time 57 is template for phase 13 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=5 [Phase 11.2] | State at t=23 [Phase 11.3] >> will be labelled as Transient(s)
#> 
#> State at time 76 is template for phase 14 
#> State at time 85 is template for phase 15 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=72 [Phase 13.3] >> will be labelled as Transient(s)
#> 
#> State at time 31 is template for phase 16 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=50 [Phase 15.1] | State at t=62 [Phase 15.2] >> will be labelled as Transient(s)
#> 
#> State at time 93 is template for phase 17 
#> ...Found state(s) already assigned to a phase in a previous iteration step:
#> ...State at t=16 [Phase 09.2] | State at t=25 [Phase 09.3] >> will be labelled as Transient(s)
#> 
#> 
#> Found 17 phases with at least 2 states.
plotRN_phaseTimeSeries(phase_out)


```
