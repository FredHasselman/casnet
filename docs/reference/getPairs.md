# Get all combinations

Get all combinations

## Usage

``` r
getPairs(x, y = NULL)
```

## Arguments

- x:

  A vector with values.

- y:

  An optional vector with values

## Value

a data frame with pairs of values

## Examples

``` r
getPairs(x=1:6)
#>    X1 X2
#> 1   1  2
#> 2   1  3
#> 3   1  4
#> 4   1  5
#> 5   1  6
#> 6   2  3
#> 7   2  4
#> 8   2  5
#> 9   2  6
#> 10  3  4
#> 11  3  5
#> 12  3  6
#> 13  4  5
#> 14  4  6
#> 15  5  6

getPairs(x=1:3, y=c("A","B","C"))
#>   Var1 Var2
#> 1    1    A
#> 2    2    A
#> 3    3    A
#> 4    1    B
#> 5    2    B
#> 6    3    B
#> 7    1    C
#> 8    2    C
#> 9    3    C
```
