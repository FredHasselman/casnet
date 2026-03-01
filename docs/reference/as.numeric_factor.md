# Numeric factor to numeric vector

Converts a factor with numeric levels to a numeric vector, using the
values of the levels.

## Usage

``` r
as.numeric_factor(x, keepNA = FALSE, sortUnique = FALSE)
```

## Arguments

- x:

  A factor based on numeric values.

- keepNA:

  Keep NA values (`TRUE`), or remove them (default = `FALSE`)

- sortUnique:

  Should the unique character values be sorted? (default = `FALSE`)

## Value

A numeric vector with factor levels as names.

## Examples

``` r
f <- factor(round(runif(10,0,9)))
as.numeric_factor(f)
#> 5 6 6 0 2 3 6 4 4 6 
#> 5 6 6 0 2 3 6 4 4 6 

# Add NAs
f <- factor(c(round(runif(9,0,9)),NA))
as.numeric_factor(f)
#> 9 2 2 6 4 6 6 1 7 
#> 9 2 2 6 4 6 6 1 7 
as.numeric_factor(f, keepNA = TRUE)
#>    9    2    2    6    4    6    6    1    7 <NA> 
#>    9    2    2    6    4    6    6    1    7   NA 

```
