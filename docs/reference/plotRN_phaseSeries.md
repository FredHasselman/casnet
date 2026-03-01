# Phase Series plot

Plot the sequence of phases as a time series.

## Usage

``` r
plotRN_phaseSeries(
  phaseOutput,
  showEpochLegend = TRUE,
  epochColours = NULL,
  epochLabel = "Phase",
  excludeOther = FALSE,
  excludeNorec = TRUE,
  excludeVars = "",
  excludePhases = "",
  plotCentroid = FALSE,
  returnGraph = FALSE
)
```

## Arguments

- phaseOutput:

  Output from function [rn_phases](rn_phases.md)

- excludeOther:

  Exclude the default Phase "Other"

- excludeNorec:

  Exclude the default Phase "No recurrence"

- excludeVars:

  Exclude specific dimension variables by name. Leave empty to include
  all variables (default = `""`)

- excludePhases:

  Exclude Phases by their name (variable `phase_name`). Leave empty to
  include all Phases (after `excludeOther` and `excludeNorec`) (default
  = `""`)

## Value

phase density plot

## Examples

``` r
RN <- rn(y1 = rnorm(100), weighted = TRUE)
#> `emRad` was set to NA due to the value `weighted = TRUE`, if you want an unthresholded matrix set `weighted = FALSE`
#> Set `weightedBy` to 'si' due to the value `weighted = TRUE`
phase_out <- rn_phases(RN)
#> This function uses an old version of the phase search algorithm...
#> Use rn_phaseInfo(), the new algorithm will produce different results!
#> 
#> Found 18 phases with at least 2 states.
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
plotRN_phaseDensity(phase_out)
#> Error in ungroup(.): could not find function "ungroup"

```
