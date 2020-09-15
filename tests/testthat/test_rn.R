context("Check rn")
library(invctr)
library(casnet)

test_that("The rn command", {
  expect_type(rn(y1 = rnorm(100)), "double")
  expect_type(rn(y1 = rnorm(100), y2 = rnorm(100)), "double")
  expect_type(rn(y1 = rnorm(100), doPlot = TRUE), "double")
})
