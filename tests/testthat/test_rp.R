context("Check rp")
library(casnet)

test_that("The rp command", {
  expect_type(rp(y1 = rnorm(100)), "S4")
  expect_type(rp(y1 = rnorm(100), y2 = rnorm(100)), "S4")
})
