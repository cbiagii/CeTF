context("PCIT")

data('simNorm')

test_that("PCIT", {
  expect_true(is.list(PCIT(simNorm)))
})
