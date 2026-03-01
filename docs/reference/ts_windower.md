# Get sliding window indices

Get sliding window indices

## Usage

``` r
ts_windower(
  y,
  win = length(y),
  step = NA,
  overlap = NA,
  adjustY = NA,
  alignment = c("r", "c", "l")[1],
  silent = TRUE
)
```

## Arguments

- y:

  A time series or numeric vector

- win:

  Size of the window to slide across `y`

- step:

  Size of steps between windows. Can be larger than `win`, but is
  ignored if `overlap` is not NA.

- overlap:

  A value between `[0 .. 1]`. If overlap is not `NA` (default), the
  value of `step` is ignored and set to `floor(overlap*win)`. This
  produces indices in which the size of `step` is always smaller than
  `win`, e.g. for fluctuation analyses that use binning procedures to
  represent time scales.

- adjustY:

  If not `NA`, or, `FALSE` a list object with fields that match one or
  more arguments of [ts_trimfill](ts_trimfill.md) (except for `x,y`),
  e.g. `list(action="trim.NA",type="end",padding=NA,silent=TRUE)`. See
  `Return value` below for details.

- alignment:

  Whether to right (`"r"`), center (`"c"`), or left (`"l"`) align the
  window.

- silent:

  Silent mode.

## Value

If `adjustY` is a list object with fields that represent arguments of
[ts_trimfill](ts_trimfill.md), then the (adjusted) vector `y` is
returned with an attribute `"windower"`. This is a list object with
fields that contain the indices for each window that fits on `y`, given
`win`, `step` or `overlap` and the settings of `adjustY`. If
`adjustY = NA`, only the list object is returned. An attribute `time` is
returned with the time index based on the `alignment` arguments.

## See also

Other Time series operations: [`ts_center()`](ts_center.md),
[`ts_changeindex()`](ts_changeindex.md),
[`ts_checkfix()`](ts_checkfix.md), [`ts_detrend()`](ts_detrend.md),
[`ts_diff()`](ts_diff.md), [`ts_discrete()`](ts_discrete.md),
[`ts_duration()`](ts_duration.md), [`ts_embed()`](ts_embed.md),
[`ts_integrate()`](ts_integrate.md), [`ts_levels()`](ts_levels.md),
[`ts_peaks()`](ts_peaks.md),
[`ts_permtest_block()`](ts_permtest_block.md),
[`ts_permtest_transmat()`](ts_permtest_transmat.md),
[`ts_rasterize()`](ts_rasterize.md), [`ts_sd()`](ts_sd.md),
[`ts_slice()`](ts_slice.md), [`ts_slopes()`](ts_slopes.md),
[`ts_standardise()`](ts_standardise.md),
[`ts_sumorder()`](ts_sumorder.md), [`ts_symbolic()`](ts_symbolic.md),
[`ts_trimfill()`](ts_trimfill.md)

Other Tools for windowed analyses: [`plotMRN_win()`](plotMRN_win.md)
