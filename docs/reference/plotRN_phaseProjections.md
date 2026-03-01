# Plot Phase Space Projection

2D umap projection of multidimensional Phase Space. Categories "No
recurrence" and "Other" are removed.

## Usage

``` r
plotRN_phaseProjections(
  RNdist,
  phaseOutput,
  epochColours = NULL,
  showEpochLegend = TRUE,
  epochLabel = "Phase"
)
```

## Arguments

- RNdist:

  A distance matrix (unthresholded) created with [rn](rn.md) or
  [rp](rp.md)

- phaseOutput:

  Output from function [rn_phaseInfo](rn_phaseInfo.md)

- PhaseOutput:

  Output from function [rn_phaseInfo](rn_phaseInfo.md) based on the
  thresholded matrix in `RNdist`

## Value

Dataframe with coordinates
