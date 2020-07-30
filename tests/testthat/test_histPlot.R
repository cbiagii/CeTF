context("histPlot")

data('simNorm')
results0 <- list()
results <- PCIT(simNorm[1:10, ])

test_that('histPlot throws an error when there are no result', {
  expect_error(histPlot(results0))
})

test_that('histPlot throws an error when the wrong element of the list is selected', {
  expect_error(histPlot(results$tab))
})
