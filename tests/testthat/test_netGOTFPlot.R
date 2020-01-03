context("netGOTFPlot")

library(org.Hs.eg.db)
data(CeTFdemo)

genes <- unique(c(as.character(NetworkData(CeTFdemo, "network1")[["gene1"]]),
                  as.character(NetworkData(CeTFdemo, "network1")[["gene2"]])))
t1 <- data.frame()
t2 <- data.frame()

test_that('netGOTFPlot throws an error when there is none netCond input', {
  expect_error(netGOTFPlot(netGO = t2,
                           keyTFs = NetworkData(CeTFdemo, "keytfs")))
})

test_that('netGOTFPlot throws an error when there is none netGO input', {
  expect_error(netGOTFPlot(netCond = NetworkData(CeTFdemo, "network1"),
                           keyTFs = NetworkData(CeTFdemo, "keytfs")))
})

test_that('netGOTFPlot throws an error when there is none keyTFs input', {
  expect_error(netGOTFPlot(netCond = NetworkData(CeTFdemo, "network1"),
                           netGO = t2))
})

