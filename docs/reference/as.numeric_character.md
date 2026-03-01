# Character vector to named numeric vector

Converts a character vector to a named numeric vector, with the
character elements as names.

## Usage

``` r
as.numeric_character(x, sortUnique = FALSE, keepNA = FALSE)
```

## Arguments

- x:

  A character vector

- sortUnique:

  Should the unique character values be sorted? (default = `FALSE`)

- keepNA:

  Keep NA values (`TRUE`), or remove them (default = `FALSE`)

## Value

A named numeric vector

## Examples

``` r
f <- letters
as.numeric_character(f)
#>  a  b  c  d  e  f  g  h  i  j  k  l  m  n  o  p  q  r  s  t  u  v  w  x  y  z 
#>  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 
```
