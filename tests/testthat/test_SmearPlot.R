context("SmearPlot")

data(CeTFdemo)

test_that('SmearPlot throws an error when there is none input for object parameter', {
  expect_error(SmearPlot(diffMethod = 'Reverter',
                           lfc = 1.5,
                           conditions = c('untrt', 'trt'),
                           TF = 'ENSG00000185917',
                           label = FALSE))
})

test_that('SmearPlot throws an error when there is none input for diffMethod parameter', {
  expect_error(SmearPlot(object = CeTFdemo,
                           lfc = 1.5,
                           conditions = c('untrt', 'trt'),
                           TF = 'ENSG00000185917',
                           label = FALSE))
})

test_that('SmearPlot throws an error when there is none input for conditions parameter', {
  expect_error(SmearPlot(object = CeTFdemo,
                           diffMethod = 'Reverter',
                           lfc = 1.5,
                           TF = 'ENSG00000185917',
                           label = FALSE))
})

test_that('SmearPlotTF throws an error when there is none input for TF parameter', {
  expect_error(SmearPlotTF(object = CeTFdemo,
                           diffMethod = 'Reverter',
                           lfc = 1.5,
                           conditions = c('untrt', 'trt'),
                           label = FALSE))
})
