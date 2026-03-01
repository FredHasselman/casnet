# Change Profile

Change Profile

## Usage

``` r
ts_cp(y, win, align = "right", keepNA = TRUE)
```

## Arguments

- y:

  Time series

- win:

  A window in which

- align:

  Alignment of the window, see
  [`zoo::rollapply()`](https://rdrr.io/pkg/zoo/man/rollapply.html)
  (default = `"right"`)

- keepNA:

  Remove `NA` or return `y` with `NA` (default = `TRUE`)

## Value

Transformed time series

## References

Hasselman, F. and Bosman, A.M.T. (2020). Studying Complex Adaptive
Systems With Internal States: A Recurrence Network Approach to the
Analysis of Multivariate Time-Series Data Representing Self-Reports of
Human Experience. Frontiers of Applied Mathematics and Statistics, 6:9.
doi: 10.3389/fams.2020.00009
