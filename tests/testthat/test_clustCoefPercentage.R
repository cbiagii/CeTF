context("clustCoefPercentage")

data('simNorm')
results <- PCIT(simNorm)

test_that("clustCoefPercentage", {
  expect_true(is.numeric(clustCoefPercentage(results$adj_sig)))
})
