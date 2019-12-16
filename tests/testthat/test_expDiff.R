context("expDiff")

data('simCounts')
anno <- data.frame(cond = c(rep('cond1', 10), rep('cond2', 10)),
                   row.names = colnames(simCounts))
colnames(simCounts) <- paste(colnames(simCounts), anno$cond, sep = "_")


test_that('expDiff throws an error when there are no anno parameter', {
  expect_error(expDiff(exp = simCounts,
                       anno = data.frame(),
                       conditions = c('cond1', 'cond2'),
                       lfc = 2,
                       padj = 0.05,
                       diffMethod = "Reverter"))
})

test_that('expDiff throws an error when there are no conditions selected', {
  expect_error(expDiff(exp = simCounts,
                       anno = anno,
                       conditions = c(),
                       lfc = 2,
                       padj = 0.05,
                       diffMethod = "Reverter"))
})

test_that("expDiff", {
  expect_true(is.matrix(simCounts) | is.data.frame(simCounts))
})

