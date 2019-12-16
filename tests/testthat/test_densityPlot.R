context("densityPlot")

data('simNorm')
results0 <- list()
results <- PCIT(simNorm)

test_that('densityPlot throws an error when there are no result', {
  expect_error(densityPlot(results0))
})

test_that('densityPlot throws an error when the wrong element of the list is selected', {
  expect_error(densityPlot(results$tab))
})
