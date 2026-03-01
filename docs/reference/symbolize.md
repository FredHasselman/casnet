# Convert numeric vectors to symbolic vectors.

`symbolize` converts numeric vectors to symbolic vectors. It is a helper
function for `muti`.

## Usage

``` r
symbolize(xy)
```

## Arguments

- xy:

  An n x 2 `matrix` or `data.frame` containing the two vectors of
  interest.

## Value

An (n-2) x 2 `matrix` of integer symbols that indicate whether the i-th
value, based on the i-1 and i+1 values, is a "trough" (=1), "decrease"
(=2), "same" (=3), "increase" (=4), or "peak" (=5).

## Author

Mark Scheuerell (https://mdscheuerell.github.io/muti/)
