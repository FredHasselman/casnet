# Recurrence Network Measures

Recurrence Network Measures

## Usage

``` r
rn_measures(g, cumulative = TRUE, silent = TRUE)
```

## Arguments

- g:

  An igraph object. If V(g)\$name is set the labels will be returned in
  a column.

- cumulative:

  Only consider out-degree.

- silent:

  Siletn(ish) mode

## Value

A list with data frames with common vertex, edge amnd global network
measures.
