context("densityPlot")

data('simNorm')
results0 <- list()
results <- PCIT(simNorm[1:10, ])

test_that('densityPlot throws an error when there are no result', {
  expect_error(densityPlot(mat1 = results0,
                          mat2 = results0,
                          threshold = 0.5))
})

test_that('densityPlot throws an error when the wrong element of the list is selected', {
  expect_error(densityPlot(mat1 = results$tab,
                          mat2 = results$adj_sig,
                          threshold = 0.5))
  
  expect_error(densityPlot(mat1 = results$adj_raw,
                          mat2 = results$tab,
                          threshold = 0.5))
  
  expect_error(densityPlot(mat1 = results$tab,
                          mat2 = results$tab,
                          threshold = 0.5))
})
