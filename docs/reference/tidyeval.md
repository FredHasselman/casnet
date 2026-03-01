# Tidy eval helpers

- [`rlang::sym()`](https://rlang.r-lib.org/reference/sym.html) creates a
  symbol from a string and
  [`rlang::syms()`](https://rlang.r-lib.org/reference/sym.html) creates
  a list of symbols from a character vector.

- [`rlang::expr()`](https://rlang.r-lib.org/reference/expr.html) and
  [`rlang::quo()`](https://rlang.r-lib.org/reference/defusing-advanced.html)
  quote one expression. `quo()` wraps the quoted expression in a
  quosure.

  The plural variants
  [`rlang::exprs()`](https://rlang.r-lib.org/reference/defusing-advanced.html)
  and
  [`rlang::quos()`](https://rlang.r-lib.org/reference/defusing-advanced.html)
  return a list of quoted expressions or quosures.

- [`rlang::enexpr()`](https://rlang.r-lib.org/reference/defusing-advanced.html)
  and [`rlang::enquo()`](https://rlang.r-lib.org/reference/enquo.html)
  capture the expression supplied as argument by the user of the current
  function (`enquo()` wraps this expression in a quosure).

  [`rlang::enexprs()`](https://rlang.r-lib.org/reference/defusing-advanced.html)
  and [`rlang::enquos()`](https://rlang.r-lib.org/reference/enquo.html)
  capture multiple expressions supplied as arguments, including `...`.

`exprs()` is not exported to avoid conflicts with `Biobase::exprs()`,
therefore one should always use
[`rlang::exprs()`](https://rlang.r-lib.org/reference/defusing-advanced.html).

To learn more about tidy eval and how to use these tools, visit
<http://rlang.r-lib.org> and the [Metaprogramming
section](https://adv-r.hadley.nz/meta.html) of [Advanced
R](https://adv-r.hadley.nz).
