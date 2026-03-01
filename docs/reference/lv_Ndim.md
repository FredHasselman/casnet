# Lotka-Volterra model for N species

The default simulates a system of `Ndim = 4` coupled competitive
Lotka-Volterra equations studied by Vano et al. (2006) using RK4
numerical integration. Vano et al. describe the dynamics of the
resulting attractor as chaotic (bounded, quasi-periodic, sensitive
dependence on initial conditions), also see [this Wiki
page](https://en.wikipedia.org/wiki/Competitive_Lotka%E2%80%93Volterra_equations#4-dimensional_example).

## Usage

``` r
lv_Ndim(
  Nsim = 1000,
  Ndim = 4,
  Y0 = rep(0.1, Ndim),
  r = c(1, 0.72, 1.53, 1.27),
  A = matrix(c(1, 1.09, 1.52, 0, 0, 1, 0.44, 1.36, 2.33, 0, 1, 0.47, 1.21, 0.51, 0.35,
    1), nrow = 4, byrow = TRUE),
  K = rep(1, Ndim),
  h = 0.5
)
```

## Arguments

- Nsim:

  How many data points to simulate (default = `1000`)

- Ndim:

  How many dimensions (coupled L-V equations) to use (default = `4`)

- Y0:

  A vector with `Ndim` elements representing the initial value (default
  = `rep(.1, Ndim)`)

- r:

  A vector with `Ndim` elements representing the growth rates for each
  dimension. The default settings are taken from Vano et al. (2006)
  (default = `c(1, 0.72, 1.53, 1.27)`)

- A:

  The interaction matrix of `Ndim X Ndim` representing the nature of the
  coupling between each dimension. `A[1,2]` will appear as a parameter
  in the equation for `Y1` setting the magnitude of the interaction with
  `Y2`. Conversely `A[2,1]` will appear in the equation for `Y2` setting
  the magnitude of the interaction with `Y1`. The diagonal elements have
  to be `1`. The default settings are taken from Vano et al. (2006)

- K:

  A vector with `Ndim` elements representing the carrying capacity for
  each dimension. (default = `rep(1, Ndim)`)

- h:

  The Euler parameter for RK4 integration.

- returnLongData:

  Return the data in long format for easy plotting (default = `FALSE`)

## Value

A matrix of `Ndim X Nsim` with simulated values.

## References

Vano, J. A., Wildenberg, J. C., Anderson, M. B., Noel, J. K., & Sprott,
J. C. (2006). Chaos in low-dimensional Lotka–Volterra models of
competition. *Nonlinearity, 19*(10), 2391.

## Examples

``` r
library(plot3D)

Y <- lv_Ndim()
#> Error in laply(1:Ndim, function(i) {    matrix(c(Y0[i], rep(NA, Nsim - 1)), nrow = 1)}): could not find function "laply"
df <- as.data.frame(apply(t(Y), 2, ts_standardise))
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 't': object 'Y' not found
lines3D(df$Y1, df$Y2, df$Y3, colvar = df$Y4, clab = "Y4", xlab= "Y1", ylab = "Y2", zlab ="Y3", ticktype = "detailed")
#> Error in df$Y4: object of type 'closure' is not subsettable
```
