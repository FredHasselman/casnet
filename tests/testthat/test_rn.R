context("Check rn")
library(invctr)
library(casnet)

test_that("The rn command", {
  expect_type(rn(y1 = rnorm(100)), "S4")
  expect_type(rn(y1 = rnorm(100), y2 = rnorm(100)), "S4")
  expect_type(rn(y1 = rnorm(100), doPlot = TRUE), "S4")
})
