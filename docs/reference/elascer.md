# Elastic Scaler - A Flexible Rescale Function

The 'elastic scaler'will rescale numeric vectors (1D, or columns in a
matrix or data.frame) to a user defined minimum and maximum, either
based on the extrema in the data, or, a minimum and maximum defined by
the user.

## Usage

``` r
elascer(
  x,
  mn = NA,
  mx = NA,
  lo = 0,
  hi = 1,
  groupwise = FALSE,
  keepNA = TRUE,
  boundaryPrecision = NA,
  tol = .Machine$double.eps^0.5
)
```

## Arguments

- x:

  Input vector or data frame.

- mn:

  Minimum value of original, defaults to `min(x, na.rm = TRUE)` if set
  to `NA`.

- mx:

  Maximum value of original, defaults to `max(x, na.rm = TRUE)` if set
  to `NA`.

- lo:

  Minimum value to rescale to, defaults to `0`.

- hi:

  Maximum value to rescale to, defaults to `1`.

- groupwise:

  If `x` is a data frame with `2+` columns, `mn = NA` and/or `mx = NA`
  and `groupwise = TRUE`, scaling will occur

- keepNA:

  Keep `NA` values?

- boundaryPrecision:

  If set to `NA` the precision of the input will be the same as the
  precision of the output. This can cause problems when detecting values
  that lie just outside of, or, exactly on boundaries given by `lo` and
  `hi`, e.g. after saving the data using a default precision. Setting
  `boundaryPrecision` to an integer value will ensure that the
  boundaries of the new scale are given by
  `round(..., digits = boundaryPrecision)`. Alternatively one could just
  round all the output after rescaling to a desired precision (default =
  `NA`)

- tol:

  The tolerance for deciding wether a value is on the boundary `lo` or
  `hi` (default = `.Machine$double.eps^0.5)`)

## Value

scaled inout

## Details

Three uses:

1.  elascer(x) - Scale x to data range: min(x.out)==0; max(x.out)==1

2.  elascer(x,mn,mx) - Scale x to arg. range: min(x.out)==mn==0;
    max(x.out)==mx==1

3.  elascer(x,mn,mx,lo,hi) - Scale x to arg. range: min(x.out)==mn==lo;
    max(x.out)==mx==hi

## Examples

``` r
# Works on numeric objects
somenumbers <- cbind(c(-5,100,sqrt(2)),c(exp(1),0,-pi))

# Using the defaults:
# 1. mn and mx are derived globally (groupWise = FALSE)
# 2. values rescaled to hi and lo are integers, 0 and 1 respectively
elascer(somenumbers)
#>           V1         V2
#> 1 0.00000000 0.07350745
#> 2 1.00000000 0.04761905
#> 3 0.06108775 0.01769912

# If the data contain values < mn they will return as < lo
elascer(somenumbers,mn=-100)
#>          V1        V2
#> 1 0.4750000 0.5135914
#> 2 1.0000000 0.5000000
#> 3 0.5070711 0.4842920
# If the data contain values > mx they will return > hi
elascer(somenumbers,mx=99)
#>           V1         V2
#> 1 0.00000000 0.07421425
#> 2 1.00961538 0.04807692
#> 3 0.06167513 0.01786930

# Effect of setting groupWise
elascer(somenumbers,lo=-1,hi=1)
#>           V1         V2
#> 1 -1.0000000 -0.8529851
#> 2  1.0000000 -0.9047619
#> 3 -0.8778245 -0.9646018
elascer(somenumbers,lo=-1,hi=1, groupwise = TRUE)
#>           V1          V2
#> 1 -1.0000000  1.00000000
#> 2  1.0000000  0.07223889
#> 3 -0.8778245 -1.00000000

elascer(somenumbers,mn=-10,mx=100,lo=-1,hi=4)
#>           V1         V2
#> 1 -0.7727273 -0.4218963
#> 2  4.0000000 -0.5454545
#> 3 -0.4811721 -0.6882542
elascer(somenumbers,mn= NA,mx=100,lo=-1,hi=4, groupwise = TRUE)
#>           V1         V2
#> 1 -1.0000000 -0.7159306
#> 2  4.0000000 -0.8477049
#> 3 -0.6945613 -1.0000000

# Effect of setting boundaryPrecision
x <- rbind(1/3, 1/7)

re1 <- elascer(x, lo = 0, hi = 1/13, boundaryPrecision = NA)
max(re1)==0.07692308 # FALSE
#> [1] FALSE
max(re1)==1/13       # TRUE
#> [1] TRUE

re2 <- elascer(x, lo = 0, hi = 1/13, boundaryPrecision = 8)
max(re2)==0.07692308 # TRUE
#> [1] TRUE
max(re2)==1/13       # FALSE
#> [1] FALSE
```
