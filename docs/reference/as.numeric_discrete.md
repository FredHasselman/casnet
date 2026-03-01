# Discrete (factor or character) to numeric vector

Converts a factor with numeric levels, or, a character vector with
numeric values to a numeric vector using
[as.numeric_factor](as.numeric_factor.md), or,
[as.numeric_character](as.numeric_character.md) respectively. If an
unnamed numeric vector is passed, it will be returned as a named numeric
vector, if this vector is continuous, it will be returned discretised
(by calling [ts_discrete](ts_discrete.md)), the labels will be rounded
to \`signif(x, digits = 4).

## Usage

``` r
as.numeric_discrete(x, keepNA = FALSE, sortUnique = FALSE)
```

## Arguments

- x:

  A factor with levels that are numeric, or, a character vector
  representing numbers.

- keepNA:

  Keep NA values (`TRUE`), or remove them (default = `FALSE`)

- sortUnique:

  Should the unique character/factor level values be sorted? (default =
  `FALSE`)

## Value

A named numeric vector with original factor levels / character values /
numeric values as names.

## Examples

``` r
# Continuous
i <- runif(10,0,9)
as.numeric_discrete(i)
#>   4.48  2.608  6.596  6.953  7.871  1.574 0.3082  2.883  3.621  1.761 
#>      3      2      5      5      5      1      1      2      3      1 

# Integer
as.numeric_discrete(round(i))
#> 4 3 7 7 8 2 0 3 4 2 
#> 4 3 7 7 8 2 0 3 4 2 

# Factor with NAs
f <- factor(c(round(runif(9,0,9)),NA))
as.numeric_discrete(f)
#> 4 1 3 9 3 6 7 2 9 
#> 4 1 3 9 3 6 7 2 9 
as.numeric_discrete(f, keepNA = FALSE)
#> 4 1 3 9 3 6 7 2 9 
#> 4 1 3 9 3 6 7 2 9 

# Character vector
c <- c("Thank","you", "for", "the flowers")
as.numeric_discrete(c)
#>       Thank         you         for the flowers 
#>           1           2           3           4 
as.numeric_discrete(c, sortUnique = TRUE)
#>       Thank         you         for the flowers 
#>           1           4           2           3 

c <- c("Thank","you", "for", "the", "flowers")
as.numeric_discrete(c)
#>   Thank     you     for     the flowers 
#>       1       2       3       4       5 
as.numeric_discrete(c, sortUnique = TRUE)
#>   Thank     you     for     the flowers 
#>       1       5       3       4       2 

```
