context("tolerance")

test_that('tolerance throws an error when the parameters a, b or are not numerics', {
  expect_error(tolerance("a", 0.5, -0.87))
})

