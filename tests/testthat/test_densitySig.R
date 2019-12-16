context("densitySig")

data('simNorm')
results0 <- list()
results <- PCIT(simNorm)

test_that('densitySig throws an error when there are no result', {
  expect_error(densitySig(mat1 = results0,
                          mat2 = results0,
                          threshold = 0.5))
})

test_that('densitySig throws an error when the wrong element of the list is selected', {
  expect_error(densitySig(mat1 = results$tab,
                          mat2 = results$adj_sig,
                          threshold = 0.5))

  expect_error(densitySig(mat1 = results$adj_raw,
                          mat2 = results$tab,
                          threshold = 0.5))

  expect_error(densitySig(mat1 = results$tab,
                          mat2 = results$tab,
                          threshold = 0.5))
})
