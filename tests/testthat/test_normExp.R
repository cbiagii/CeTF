context("normExp")

data('simCounts')
tpm <- apply(simCounts, 2, function(x) {
  (1e+06 * x)/sum(x)
})


test_that("normExp", {
  expect_true(is.matrix(normExp(tpm)))
})
