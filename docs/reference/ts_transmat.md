# Transition matrix

Create a transition matrix from a discrete time series, e.g. to generate
Monte Carlo simulations.

## Usage

``` r
ts_transmat(yd, nbins = length(unique(yd)))
```

## Arguments

- yd:

  A discrete numeric vector or time series, e.g. transformed using
  [`ts_discrete()`](ts_discrete.md), or,
  [`ts_symbolic()`](ts_symbolic.md).

- nbins:

  The number of bins used to transform a continuous time series, or, the
  number of expected (given `nbins`, or, theoretically possible) values
  for a discrete series (default = `length(unique(yd))`)

## Value

A transition probability matrix

## Examples

``` r
set.seed(4321)

# Random uniform numbers
y  <- runif(10,0,20)

# Discrete version
yd <- ts_discrete(y, nbins = 10)

# Transition probabilities
ts_transmat(yd = yd, nbins = 10)
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]      [,9]     [,10]
#>  [1,]  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0 1.0000000 0.0000000
#>  [2,]  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1 0.1000000 0.1000000
#>  [3,]  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1 0.1000000 0.1000000
#>  [4,]  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0 0.0000000 1.0000000
#>  [5,]  0.5  0.0  0.0  0.0  0.0  0.5  0.0  0.0 0.0000000 0.0000000
#>  [6,]  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1 0.1000000 0.1000000
#>  [7,]  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1 0.1000000 0.1000000
#>  [8,]  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1 0.1000000 0.1000000
#>  [9,]  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0 0.6666667 0.3333333
#> [10,]  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0 0.0000000 0.0000000

# Note: The number of 'observed' bins differs from 'expected' bins
table(yd)
#> yd
#>  1  4  5  6  9 10 
#>  1  1  2  1  3  2 

# Not specifying the expected bins will generate a different matrix!
ts_transmat(yd = yd, nbins = length(unique(yd)))
#>           [,1]      [,2]      [,3]      [,4]      [,5]      [,6]
#> [1,] 0.1666667 0.1666667 0.1666667 0.1666667 0.1666667 0.1666667
#> [2,] 0.1666667 0.1666667 0.1666667 0.1666667 0.1666667 0.1666667
#> [3,] 0.1666667 0.1666667 0.1666667 0.1666667 0.1666667 0.1666667
#> [4,] 0.1666667 0.1666667 0.1666667 0.1666667 0.1666667 0.1666667
#> [5,] 0.5000000 0.0000000 0.0000000 0.0000000 0.0000000 0.5000000
#> [6,] 0.1666667 0.1666667 0.1666667 0.1666667 0.1666667 0.1666667
```
