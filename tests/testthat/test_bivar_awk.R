context("bivar.awk")

df <- data.frame(a = sample(1:1000, 100, replace=TRUE),
                 b = sample(1:1000, 100, replace=TRUE))

test_that("bivar.awk", {
  expect_true('mean1' %in% names(bivar.awk(df)))
  expect_true('mean2' %in% names(bivar.awk(df)))
})
