context("Check crqa_cl")
library(invctr)
library(casnet)

test_that("Commandline crp", {
  expect_output(str(crqa_cl(y1 = rnorm(100))), "Performing auto-RQA", fixed = TRUE)
  expect_output(str(crqa_cl(y1 = rnorm(100), y2 = rnorm(100))), "Performing cross-RQA", fixed = TRUE)
  expect_output(str(crqa_cl(y1 = rnorm(100), win=50)), "Calculating recurrence measures in a window", fixed = TRUE)
})
