# Create Cauchy Flight

Creates a Cauchy flight by taking increments from the Cauchy
distributions implemented as the stable distribution
([`stabledist::rstable()`](https://rdrr.io/pkg/stabledist/man/dist-stable.html))
with index paramter `alpha = 1` and skewness parameter `beta = 0`.

## Usage

``` r
flight_Cauchy(
  N = 1000,
  ndims = 2,
  alpha = 1,
  beta = 0,
  scale = 1,
  location = 0
)
```

## Arguments

- N:

  Length of time series (default = `1000`)

- ndims:

  Number of dimensions (default = `2`)

- alpha:

  Index of stability parameter in `(0,2]`

- beta:

  Skewness parameter in `[-1,1]`

- scale:

  Scale parameterin `(0,Inf)`

- location:

  Location (shift) parameter in `[-Inf,Inf]`

## Value

A data frame with `ndims` columns and `N` rows.

## Examples

``` r
df <- flight_Cauchy()
plot(density(diff(df$dim1)))

plot(df$dim1, df$dim2, type = "l")

```
