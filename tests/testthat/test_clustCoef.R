context("clustCoef")

data('simNorm')
results <- PCIT(simNorm)

test_that("clustCoef", {
  expect_true(is.numeric(clustCoef(results$adj_sig)))
})
