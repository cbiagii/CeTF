context("countsToTPM")

data('simCounts')
tpm <- countsToTPM(simCounts)

test_that("countsToTPM", {
  expect_true(is.matrix(tpm))
})
