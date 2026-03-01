# Plot Phase Space Projection

2D umap projection of multidimensional Phase Space. Categories "No
recurrence" and "Other" are removed.

## Usage

``` r
plotRN_phaseProjection(
  RNdist,
  phaseOutput,
  epochColours = NULL,
  showEpochLegend = TRUE,
  epochLabel = "Phase",
  excludeOther = TRUE,
  excludeNorec = TRUE
)
```

## Arguments

- RNdist:

  A distance matrix (unthresholded) created with [rn](rn.md) or
  [rp](rp.md)

- phaseOutput:

  Output from function [rn_phases](rn_phases.md)

- excludeOther:

  Exclude the default Phase "Other"

- excludeNorec:

  Exclude the default Phase "No recurrence"

- PhaseOut:

  Output from function [rn_phases](rn_phases.md) based on the
  thresholded matrix in `RNdist`

## Value

Dataframe with coordinates
