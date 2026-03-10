# Prepare time series for fluctuation analysis

Prepare time series for fluctuation analysis

## Usage

``` r
fd_prepSeries(
  y,
  fs = NULL,
  removeTrend = c("no", "poly", "adaptive", "bridge")[2],
  polyOrder = 1,
  standardise = c("none", "mean.sd", "median.mad")[2],
  adjustSumOrder = FALSE,
  scaleMin = NA,
  scaleMax = NA,
  scaleResolution = NA,
  scaleS = NA,
  Nyquist = TRUE,
  overlap = NA,
  silent = TRUE
)
```

## Arguments

- y:

  y

- fs:

  fs

- removeTrend:

  rt

- polyOrder:

  po

- standardise:

  st

- adjustSumOrder:

  ao

- scaleMin:

  sm

- scaleMax:

  smm

- scaleResolution:

  sr

- scaleS:

  sc

- Nyquist:

  Ny

- overlap:

  ov

- silent:

  si

## Value

transformed series
