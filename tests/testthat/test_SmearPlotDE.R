context("SmearPlotDE")

data(pcitrifExample)

test_that('SmearPlotDE throws an error when there is none input for object parameter', {
  expect_error(SmearPlotDE(diffMethod = 'Reverter',
                           lfc = 1.5,
                           conditions = c('untrt', 'trt')))
})

test_that('SmearPlotDE throws an error when there is none input for diffMethod parameter', {
  expect_error(SmearPlotDE(object = pcitrifExample,
                           lfc = 1.5,
                           conditions = c('untrt', 'trt')))
})

test_that('SmearPlotDE throws an error when there is none input for conditions parameter', {
  expect_error(SmearPlotDE(object = pcitrifExample,
                           diffMethod = 'Reverter',
                           lfc = 1.5))
})
