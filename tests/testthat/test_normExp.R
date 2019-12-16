context("normExp")

data('simCounts')
tpm <- countsToTPM(simCounts)

test_that("normExp", {
  expect_true(is.matrix(normExp(tpm)))
})
