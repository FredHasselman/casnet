# Prepare matrix

Prepare matrix

## Usage

``` r
rp_prep(
  RP,
  emRad = NULL,
  DLmin = 2,
  VLmin = 2,
  HLmin = 2,
  DLmax = length(Matrix::diag(RP)),
  VLmax = length(Matrix::diag(RP)),
  HLmax = length(Matrix::diag(RP)),
  AUTO = FALSE,
  chromatic = FALSE,
  matrices = FALSE,
  doHalf = FALSE
)
```

## Arguments

- RP:

  Recurrence plot

- emRad:

  Radiuc

- DLmin:

  Minimal diagonal line length

- VLmin:

  Minimal vertical line length

- HLmin:

  Minimal horizontal line length

- DLmax:

  Maximal diagonal line length

- VLmax:

  Maximal vertical line length

- HLmax:

  Maximal horizontal line length

- AUTO:

  Is this an AUTO RQA?

- chromatic:

  Chromatic RQA?

- matrices:

  Return Matrices?

- doHalf:

  Analyse half of the matrix?

## Value

A prepped matrix
