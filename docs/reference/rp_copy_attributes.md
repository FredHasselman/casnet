# Copy Matrix Attributes

Simple attribute copy used in `casnet` to convert between `matrix` and
`Matrix` classes and back.

## Usage

``` r
rp_copy_attributes(
  source,
  target,
  source_remove = c("names", "row.names", "class", "dim", "dimnames", "x")
)
```

## Arguments

- source:

  Source matrix

- target:

  Target matrix

- source_remove:

  Remove these attribute fields from the source before copying.

## Value

The target matrix with attributes copied from the source matrix.
