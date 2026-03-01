# Fast (C)RQA (command line crp)

This function will run the [commandline Recurrence
Plots](http://tocsy.pik-potsdam.de/commandline-rp.php) executable
provided by Norbert Marwan.

## Usage

``` r
rp_cl(
  y1,
  y2 = NULL,
  emDim = 1,
  emLag = 1,
  emRad = NA,
  DLmin = 2,
  VLmin = 2,
  theiler = 0,
  win = min(length(y1), ifelse(is.null(y2), (length(y1) + 1), length(y2)), na.rm = TRUE),
  step = win,
  JRP = FALSE,
  distNorm = c("EUCLIDEAN", "MAX", "MIN", "OP")[[1]],
  standardise = c("none", "mean.sd", "median.mad")[1],
  returnMeasures = TRUE,
  returnRPvector = FALSE,
  returnLineDist = FALSE,
  doPlot = c("noplot", "rp", "distmat")[[1]],
  path_to_rp = getOption("casnet.path_to_rp"),
  saveOut = FALSE,
  path_out = NULL,
  file_ID = NULL,
  silent = TRUE,
  surrogateTest = FALSE,
  targetValue = 0.05,
  useParallel = FALSE,
  ...
)
```

## Arguments

- y1:

  Time series 1

- y2:

  Time series 2 for Cross Recurrence Analysis (default = `NULL`)

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

- standardise:

  Standardise data: `"none"` (default), `"mean.sd"`, or `"median.mad"`

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

- surrogateTest:

  Perform surrogate tests. If `TRUE`, will run surrogate tests using
  default settings for a two-sided test of \\H_0: The data generating
  process is a rescaled linear Gaussian process\\ at \\\alpha = .05\\
  (arguments `ns = 39, fft = TRUE, amplitude = TRUE`)

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

## Value

A list object containing 1-3 elements, depending on arguments requesting
output.

- `rqa_measures` - A list of the (C)RQA measures returned if
  `returnMeasures = TRUE`:

  1.  RR = Recurrence rate

  2.  DET = Determinism

  3.  DET_RR = Ratio DET/RR

  4.  LAM = Laminarity

  5.  LAM_DET = Ratio LAM/DET

  6.  L_max = maximal diagonal line length

  7.  L_mean = mean diagonal line length

  8.  L_entr = Entropy of diagonal line length distribution

  9.  DIV = Divergence (1/L_max)

  10. V_max = maximal vertical line length

  11. TT = Trapping time

  12. V_entr = Entropy of vertical line length distribution

  13. T1 = Recurrence times 1st type

  14. T2 = Recurrence times 2nd type

  15. W_max = Max interval length

  16. W_mean = Mean of interval lengths

  17. W_entr = Entropy of interval length distribution

  18. W_prob = Probability of interval

  19. F_min = F min

- `rqa_rpvector` - The radius thresholded distance matrix (recurrence
  matrix), which can be visualised as a recurrence plot by calling
  [`rp_plot()`](rp_plot.md). If a sliding window analysis is conducted
  this will be a list of matrices and could potentially grow too large
  to handle. It is recommended you save the output to disk by setting
  `saveOut = TRUE`.

- `rqa_diagdist` - The distribution of diagonal line lengths

## Details

The `rp` executable is installed when the function is called for the
first time and is renamed to `rp`, from a platform specific filename
downloaded from http://tocsy.pik-potsdam.de/commandline-rp.php or
extracted from an archive located in the directory:
`...\\casnet\\commandline_rp\\` The file is copied to the directory:
`...\\casnet\\exec\\` The latter location is stored as an option and can
be read by calling `getOption("casnet.path_to_rp")`.

## Note

The platform specific `rp` command line executables were created by
Norbert Marwan and obtained under a Creative Commons License from the
website of the Potsdam Institute for Climate Impact Research at
<http://tocsy.pik-potsdam.de/>.

The full copyright statement on the website is as follows:

\(C\) 2004-2017 SOME RIGHTS RESERVED

University of Potsdam, Interdisciplinary Center for Dynamics of Complex
Systems, Germany

Potsdam Institute for Climate Impact Research, Transdisciplinary
Concepts and Methods, Germany

This work is licensed under a [Creative Commons
Attribution-NonCommercial-NoDerivs 2.0 Germany
License](https://creativecommons.org/licenses/by-nc-nd/2.0/de/).

More information about recurrence analysis can be found on the
[Recurrence Plot](http://www.recurrence-plot.tk) website.

## Troubleshooting

Some notes on resolving errors with `rp`.The script will first try to
download the correct executable, if that fails, the copy will have...
failed. It should be relatively easy to get `rp_cl()` working though, by
using some custom settings:

- *Copy failed* - Every time the function `rp_cl()` is called it will
  check whether a log file `rp_instal_log.txt` is present in the
  `...\\casnet\\exec\\` directory. If you delete the `rp_instal_log.txt`
  file, and call the function, another attempt will be made to download
  a copy of the executable.

- *Copy still fails and/or no permission to copy* - If you cannot acces
  the directory `...\\casnet\\commandline_rp\\`, download the
  appropriate executable from the [Commandline Recurrence
  Plots](http://tocsy.pik-potsdam.de/commandline-rp.php) page and copy
  to a directory you do have the rights to: *execute* commands, *write*
  and *read* files. Make sure you rename the file to `rp` (`rp.exe` on
  Windows OS). Then, either pass the path to `rp` as the argument
  `path_to_rp` in the `rp_cl(.., path_to_rp = "YOUR_PATH")` function
  call, or, as a more permanent solution, set the `path_to_rp` option by
  calling `options(casnet.path_to_rp="YOUR_PATH")`.

- *Error in execution of `rp`* - This can have a variety of causes, the
  `rp` executable is called using
  [`system2()`](https://rdrr.io/r/base/system2.html) and makes use of
  the [`normalizePath()`](https://rdrr.io/r/base/normalizePath.html)
  function with argument `mustWork = FALSE`. Problems caused by specific
  OS, machine, or, locale problems (e.g. the `winslash` can be reported
  as an [issue on
  Github](https://github.com/FredHasselman/casnet/issues)). One
  execution error that occurs when the OS is not recognised properly can
  be resolved by chekcing `getOption("casnet.rp_prefix")`. On Windows OS
  this should return an empty character vector, on Linux or macOS it
  should return `"./"`. You can manually set the correct prefix by
  calling `options(casnet.rp_prefix="CORRECT OS PREFIX")` and fill in
  the prefix that is correct for your OS

## See also

Other Recurrence Quantification Analysis:
[`rp_measures()`](rp_measures.md),
[`rp_measures_main()`](rp_measures_main.md)
