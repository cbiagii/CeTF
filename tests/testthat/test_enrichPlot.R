context("enrichPlot")

data(CeTFdemo)
genes <- unique(c(as.character(NetworkData(CeTFdemo, "network1")[["gene1"]]),
                  as.character(NetworkData(CeTFdemo, "network1")[["gene2"]])))
 

test_that('enrichPlot throws an error when there is none res argument', {
  expect_error(enrichPlot(showCategory = 10, 
                          type = "circle"))
})

test_that('enrichPlot throws an error when showCategory argument ir equal 1', {
  expect_error(enrichPlot(res = data.frame(),
                          showCategory = 1, 
                          type = "circle"))
})
