context("singleNetworkPlot")

data('simCounts')
out <- NULL

test_that('singleNetworkPlot throws an error when the input is different from a pcitrif class object', {
  expect_error(singleNetworkPlot(out))
})
