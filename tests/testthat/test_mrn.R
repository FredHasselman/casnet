context("Check mrn")
library(igraph)
library(casnet)

test_that("The mrn command", {
  expect_type(mrn(layers = list(make_ring(10),make_ring(10))), "list")
  expect_type(mrn(list(make_ring(10),make_ring(10)),win = 5), "list")
})
