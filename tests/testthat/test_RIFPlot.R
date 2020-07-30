context("RIFPlot")

data(CeTFdemo)

test_that('RIFPlot throws an error when there is none input for object parameter', {
  expect_error(RIFPlot(color  = 'darkblue', 
                       type   = 'DE'))
})
