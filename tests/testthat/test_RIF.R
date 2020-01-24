context("RIF")

data('RIF_input')

test_that('RIF throws an error when there are no nta or ntf parameters', {
  expect_error(RIF(input = RIF_input,
                   nta = NULL,
                   ntf = NULL,
                   ncond1 = 10,
                   ncond2 = 10))
})

test_that('RIF throws an error when there are no ncond1 or ncond2 parameters', {
  expect_error(RIF(input = RIF_input,
                   nta = 104,
                   ntf = 50,
                   ncond1 = NULL,
                   ncond2 = NULL))
})

test_that("RIF", {
  expect_true(is.matrix(RIF_input) | is.data.frame(RIF_input))
})

