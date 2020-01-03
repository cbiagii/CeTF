context("clustCoef")

data('simNorm')
results <- PCIT(simNorm[1:10, ])

test_that("clustCoef", {
  expect_true(is.numeric(clustCoef(results$adj_sig)))
})
