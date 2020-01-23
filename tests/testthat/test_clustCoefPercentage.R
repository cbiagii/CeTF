context("clustCoefPercentage")

data('simNorm')
results <- PCIT(simNorm[1:10, ])

test_that("clustCoefPercentage", {
  expect_true(is.numeric(clustCoefPercentage(results$adj_sig)))
})
