# Diagonal Recurrence Profile

Diagonal Recurrence Profile

## Usage

``` r
rp_diagProfile(
  RM,
  diagWin = NULL,
  xname = "X-axis",
  yname = "Y-axis",
  theiler = 0,
  DLmin = 2,
  VLmin = 2,
  HLmin = 2,
  DLmax = length(Matrix::diag(RM)),
  VLmax = length(Matrix::diag(RM)),
  HLmax = length(Matrix::diag(RM)),
  doShuffle = FALSE,
  y1 = NULL,
  y2 = NULL,
  shuffleWhich = "y1",
  Nshuffle = 19,
  doEmbed = TRUE,
  AUTO = NULL,
  chromatic = FALSE,
  matrices = FALSE,
  doPlot = TRUE,
  minY = 0.5,
  returnOnlyPlot = FALSE
)
```

## Arguments

- RM:

  A binary recurrence matrix

- diagWin:

  Window around the line of synchrony

- xname:

  Label for x-axis

- yname:

  Label for y-axis

- theiler:

  Use a `theiler` window around the main diagonal (Line of
  Identity/Synchronisation) to remove auto-correlations at short
  time-lags:

  - `0` will include the main diagonal in all RQA measure calculations.

  - `1` will remove the main diagonal from all RQA measure calculations.

  - `NA` (default), will check if the matrix is symmetrical , if so, it
    will remove the diagonal by setting `theiler = 1` (Line of Identity,
    Auto-RQA), if it is not symmetrical (Line of Synchronisation,
    Cross-RQA) it will set `theiler = 0`.

  - A value greater than `1` will remove that many diagonals around and
    including the diagonal from all RQA measure calculations. So
    `theiler = 2` means exclude `2` diagonals around the main diagonal,
    including the main diagonal itself: `[-1,0,1]`. If `theiler` is a
    numeric vector of `length(theiler) == 2` it is possible to exclude
    an asymmetrical window. The values are interpreted as end points in
    a sequence of diagonal ID's, e.g. `theiler = c(-1,5)` will exclude
    `[-1,0,1,2,3,4,5]`. If `length(theiler) > 2`, the values will be
    considered individual diagonal ID's, e.g.
    `theiler = c(-3,-1,0,2,5)`, will exclude only those specific ID's.
    Also see the note.

- DLmin:

  Minimal diagonal line length (default = `2`)

- VLmin:

  Minimal vertical line length (default = `2`)

- HLmin:

  Minimal horizontal line length (default = `2`)

- DLmax:

  Maximal diagonal line length (default = length of diagonal -1)

- VLmax:

  Maximal vertical line length (default = length of diagonal -1)

- HLmax:

  Maximal horizontal line length (default = length of diagonal -1)

- doShuffle:

  Should a shuffled baseline be calculated (default = `FALSE`)

- y1:

  The original `y1` time series

- y2:

  The original `y2` time series

- shuffleWhich:

  Which of the time series should be shuffled: 'y1' or 'y2'? (default =
  'y2')

- Nshuffle:

  How many shuffled versions to make up the baseline? The default is
  `19`, which is the minimum for a one-sided surrogate test.

- doEmbed:

  If `doShuffle = TRUE`, should the data in y1 and y2 be considered
  embedded time series? The temporal order of all columns in `y2` will
  be randomly shuffled in the same way, keeping coordinates together
  (default = `FALSE`)

- AUTO:

  Auto-recurrence? (default = `FALSE`)

- chromatic:

  Force chromatic RQA? (default = `FALSE`)

- matrices:

  Return matrices? (default = `FALSE`)

- doPlot:

  Plot (default = `TRUE`)

- minY:

  The minimum value of the Y-axis. The upper Y limit is determined by
  `max(minY,max(RR))`. Set to 1 to always show the theoretical maximum
  (default = `0.5`)

- returnOnlyPlot:

  Don't plot to graphics device, but do return the plot (default =
  `FALSE`)

## Value

A plot and/or the data for the plot
