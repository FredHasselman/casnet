# mono Hurst

mono Hurst

## Usage

``` r
monoH(
  TSm,
  scaleS,
  removeTrend = c("no", "poly", "adaptive", "bridge")[2],
  polyOrder = 1,
  overlap = overlap,
  returnPLAW = FALSE,
  returnSegments = FALSE,
  removeRMSbelow = .Machine$double.eps
)
```

## Arguments

- TSm:

  TS matrix with 2 columns `t` (1st) and `y` (second)

- scaleS:

  If not `NA`, it should be a numeric vector listing the scales on which
  to evaluate the detrended fluctuations. Arguments
  `scaleMax, scaleMin, scaleResolution` will be ignored (default = `NA`)

- removeTrend:

  Method to use for global detrending (default = `"poly"`)

- polyOrder:

  Order of global polynomial trend to remove if `removeTrend = "poly"`.
  If `removeTrend = "adaptive"` polynomials `1` to `polyOrder` will be
  evaluated and the best fitting curve (R squared) will be removed
  (default = `1`)

- overlap:

  A number in `[0 ... 1]` representing the amount of 'bin overlap' when
  calculating the fluctuation. This reduces impact of arbitrary time
  series begin and end points. If `length(y) = 1024` and overlap is
  `.5`, a scale of `4` will be considered a sliding window of size `4`
  with step-size `floor(.5 * 4) = 2`, so for scale `128` step-size will
  be `64` (default = `NA`)

- returnPLAW:

  Return the power law data (default = `FALSE`)
