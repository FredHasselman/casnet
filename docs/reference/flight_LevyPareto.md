# Create a Levy-Pareto flight

Creates a Rayleigh flight by taking increments from the Normal
distributions implemented as the stable distribution
([`stabledist::rstable()`](https://rdrr.io/pkg/stabledist/man/dist-stable.html))
with index paramter `alpha = 1.5` and skewness parameter `beta = 0`.

## Usage

``` r
flight_LevyPareto(
  N = 1000,
  ndims = 2,
  alpha = 1.5,
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

## Details

Note that the increments are not strictly from the distribution called
**the** Levy distribution, but rather **a** a Levy-with-Pareto-tail-type
distribution (i.e. when `1 < alpha < 2`). Use `alpha = 1/2` and
`beta = 1` if **the** Levy distribution is required.

## Examples

``` r
# Levy-Pareto
df <- flight_LevyPareto()
plot(density(diff(df$dim1)))

plot(df$dim1, df$dim2, type = "l")


# "The" Levy distribution
df <- flight_LevyPareto(alpha = 1/2, beta = 1)
#> Warning: Not a Levy-Pareto distribution if alpha <= 1 or alpha >= 2
plot(density(diff(df$dim1)))

plot(df$dim1, df$dim2, type = "l")

```
