context("runAnalysis")

data('simCounts')

test_that('runAnalysis throws an error when there are no conditions selected', {
  expect_error(runAnalysis(mat = simCounts,
                           conditions=c(),
                           lfc = 2.57,
                           padj = 0.05,
                           TFs = paste0('TF_', 1:1000),
                           ncond1 = 10,
                           ncond2= 10,
                           tolType = 'mean',
                           diffMethod = 'Reverter'))
})

test_that('runAnalysis throws an error when there are no TFs inputed', {
  expect_error(runAnalysis(mat = simCounts,
                           conditions=c('cond1', 'cond2'),
                           lfc = 2.57,
                           padj = 0.05,
                           TFs = c(),
                           ncond1 = 10,
                           ncond2= 10,
                           tolType = 'mean',
                           diffMethod = 'Reverter'))
})

test_that('runAnalysis throws an error when there are no ncond1 or ncond2 parameters defined', {
  expect_error(runAnalysis(mat = simCounts,
                           conditions=c('cond1', 'cond2'),
                           lfc = 2.57,
                           padj = 0.05,
                           TFs = paste0('TF_', 1:1000),
                           ncond1 = NULL,
                           ncond2= NULL,
                           tolType = 'mean',
                           diffMethod = 'Reverter'))
})
