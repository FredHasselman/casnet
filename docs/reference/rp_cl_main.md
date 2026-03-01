# rp_cl_main

rp_cl_main

## Usage

``` r
rp_cl_main(
  data,
  emDim = 1,
  emLag = 1,
  emRad = NA,
  DLmin = 2,
  VLmin = 2,
  theiler = 0,
  win = min(length(y1), ifelse(is.null(y2), (length(y1) + 1), length(y2)), na.rm = TRUE),
  step = step,
  JRP = FALSE,
  distNorm = c("EUCLIDEAN", "MAX", "MIN", "OP")[[1]],
  returnMeasures = TRUE,
  returnRPvector = FALSE,
  returnLineDist = FALSE,
  doPlot = c("noplot", "rp", "distmat")[[1]],
  path_to_rp = getOption("casnet.path_to_rp"),
  saveOut = FALSE,
  path_out = NULL,
  file_ID = NULL,
  silent = TRUE,
  targetValue = 0.05,
  useParallel = FALSE,
  ...
)
```

## Arguments

- emDim:

  Embedding dimensions (default = `1`)

- emLag:

  Embedding lag (default = `1`)

- emRad:

  Radius on distance matrix (default = `1`)

- DLmin:

  Minimum length of diagonal structure to be considered a line (default
  = `2`)

- VLmin:

  Minimum length of vertical structure to be considered a line (default
  = `2`)

- theiler:

  Theiler window (default = `0`)

- win:

  Window to calculate the (C)RQA (default = minimum of length of `y1` or
  `y2`)

- step:

  Stepsize for sliding windows (default = size of `win`, so no sliding
  window)

- JRP:

  Wether to calculate a Joint Recurrence Plot (default = `FALSE`)

- distNorm:

  One of "EUCLIDEAN" (default), `"MAX", "MIN"`, or `"OP"` for an Order
  Pattern recurrence matrix

- returnMeasures:

  Return the (C)RQA measures? (default = `TRUE`)

- returnRPvector:

  Return the recurrent points in a dataframe? (default = `FALSE`)

- returnLineDist:

  Return the distribution of diagonal and horizontal line length
  distances (default = `FALSE`)

- doPlot:

  Produce a plot of the recurrence matrix by calling
  [`rp_plot()`](rp_plot.md), values can be `"rp"` (the thresholded
  recurrence matrix),`"distmat"` (the unthresholded recurrence matrix)
  or `"noplot"` (default = `"noplot"`)

- path_to_rp:

  Path to the command line executable (default = path set during
  installation, use `getOption("casnet.path_to_rp")` to see)

- saveOut:

  Save the output to files? If `TRUE` and `path_out = NA`, the current
  working directory will be used (default = `FALSE`)

- path_out:

  Path to save output if `saveOut = TRUE` (default = `NULL`)

- file_ID:

  A file ID which will be a prefix to to the filename if
  `saveOut = TRUE` (default = `NULL`, an integer will be added tot the
  file name to ensure unique files)

- silent:

  Do not display any messages (default = `TRUE`)

- targetValue:

  A value passed to `est_radius(...,type="fixed", targetMeasure="RR")`
  if `is.na(emRad)==TRUE`. This is useful for windowed analysis, it will
  estimate a new radius for each window.

- useParallel:

  Speed up calculations by using the parallel processing options
  provided by `parallel` to assign a seperate process/core for each
  window in windowed (C)RQA analysis and
  [`parallel::detectCores()`](https://rdrr.io/r/parallel/detectCores.html)
  with`logical = TRUE` to decide on the available cores (default =
  `FALSE`)

- ...:

  Additional parameters (currently not used)
