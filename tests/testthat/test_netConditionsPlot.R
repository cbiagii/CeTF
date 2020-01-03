context("netConditionsPlot")

data('simCounts')
out <- NULL

test_that('netConditionsPlot throws an error when the input is different from a CeTF class object', {
  expect_error(netConditionsPlot(out))
})
