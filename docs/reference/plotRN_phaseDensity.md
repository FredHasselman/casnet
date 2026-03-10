# Phase Density for each dimension

Create a ridge plot (requires package
[ggridges::ggridges](https://wilkelab.org/ggridges/reference/ggridges-package.html))

## Usage

``` r
plotRN_phaseDensity(
  phaseOutput,
  plotCentroid = NA,
  excludeOther = FALSE,
  excludeNorec = FALSE,
  excludeVars = "",
  excludePhases = "",
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

data frame with phase series.

## Examples

``` r
RN <- rn(cumsum(rnorm(100)), emDim = 1, emLag = 1, emRad = NA, weighted = TRUE)
#> `emRad` was set to NA due to the value `weighted = TRUE`, if you want an unthresholded matrix set `weighted = FALSE`
#> Set `weightedBy` to 'si' due to the value `weighted = TRUE`
outPhases <- rn_phases(RN)
#> This function uses an old version of the phase search algorithm...
#> Use rn_phaseInfo(), the new algorithm will produce different results!
#> 
#> Found 20 phases with at least 2 states.
plotRN_phaseDensity(outPhases)
#> Error in ungroup(.): could not find function "ungroup"
```
