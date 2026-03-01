# Course grain a matrix for plotting

Course grain a matrix for plotting

## Usage

``` r
mat_coursegrain(
  RM,
  target_height = round(NROW(RM)/2),
  target_width = round(NCOL(RM)/2),
  summary_func = NA,
  recurrence_threshold = NA,
  categorical = NA,
  output_type = 0,
  n_core = 1,
  silent = FALSE
)
```

## Arguments

- RM:

  A (recurrence) matrix

- target_height:

  How many rows? (default = `NROW(RM)/2`)

- target_width:

  How many columns? (default = `NCOL(RM)/2`)

- summary_func:

  How to summarise the values in subset `X` of `RM`. If set to `NA`, the
  function will try to pick a summary function based on the cell values:
  If `RM` is a distance matrix, `mean(X, na.rm = TRUE)` will be used; If
  it is a binary matrix
  `ifelse(mean(X, na.rm = TRUE)>recurrence_threshold,1,0)`, a
  categorical matrix (`categorical = TRUE`, or, matrix attribute
  `chromatic = TRUE`) will pick the most frequent category in the subset
  `attributes(ftable(X))$col.vars$x[[which.max(ftable(X))]]`. (default =
  `NA`)

- recurrence_threshold:

  For a binary matrix the mean of the cells to be summarised will vary
  between `0` and `1`, which essentially represents the recurrence rate
  for that subset of the matrix. If `NA` the threshold will be set to a
  value that in most cases should return a plot with a similar `RR` as
  the original plot. (default = `NA`)

- categorical:

  If set to `TRUE`, will force `summary_func` to select the most
  frequent value. If `NA` the matrix attribute `chromatic` will be used.
  If `chromatic` is not present, all values in the matrix have to be
  whole numbers as determined by
  [`plyr::is.discrete()`](https://rdrr.io/pkg/plyr/man/is.discrete.html).
  (default = `NA`)

- output_type:

  The output format for `plyr::vapply()`. (default = `0.0`)

- n_core:

  Number of cores for parallel processing. Set to `NA` to automatically
  choose cores. (default = `1`)

- silent:

  Silt-ish mode (default = `FALSE`)

## Value

A coursegrained matrix of size `target_width` by `target_height`.

## Note

This code was inspired by code published in a blog post by Guillaume
Devailly on 29-04-2020
(https://gdevailly.netlify.app/post/plotting-big-matrices-in-r/)

## Examples

``` r
# Continuous
RMc1 <- rp(cumsum(rnorm(200)))
rp_plot(RMc1)
#> Warning: `aes_()` was deprecated in ggplot2 3.0.0.
#> ℹ Please use tidy evaluation idioms with `aes()`
#> ℹ The deprecated feature was likely used in the casnet package.
#>   Please report the issue at <https://github.com/FredHasselman/casnet/issues>.

RMc2 <- mat_coursegrain(RMc1)
#> Continuous matrix... using summary function 'mean(x, na.rm = TRUE)' for coursegraining.
rp_plot(RMc2)


# Binary
RMb1 <- rp(cumsum(rnorm(200)), emRad = NA)
rp_plot(RMb1, plotMeasures = TRUE)
#> Auto-RQA, not including diagonal, theiler set to 1...

# Reported RQA measures in rp_plot will be based on the full matrix
rp_plot(RMb1, maxSize = 100^2, plotMeasures = TRUE)
#> Auto-RQA, not including diagonal, theiler set to 1...
#> NOTE: To speed up the plotting process, the RP will represent a coursegrained matrix. Set argument 'courseGrain = FALSE' to see the full matrix.
#> Auto-RQA, not including diagonal, theiler set to 1...
#> Binary matrix... using summary function 'ifelse(mean(x, na.rm = TRUE)>recurrence_threshold,1,0)' for coursegraining.

# Plotting the coursegrained matrix itself will yield different values
RMb2 <- mat_coursegrain(RMb1)
#> Auto-RQA, not including diagonal, theiler set to 1...
#> Binary matrix... using summary function 'ifelse(mean(x, na.rm = TRUE)>recurrence_threshold,1,0)' for coursegraining.
rp_plot(RMb2, plotMeasures = TRUE)
#> Auto-RQA, not including diagonal, theiler set to 1...


# Categorical
RMl1 <- rp(y1 = round(runif(100, min = 1, max = 3)), chromatic = TRUE)
rp_plot(RMl1)

RMl2 <- mat_coursegrain(RMl1, categorical = TRUE)
#> Categorical matrix... using summary function 'attributes(ftable(x))$col.vars[[which.max(ftable(x))]]' for coursegraining.
rp_plot(RMl2)

```
