context("heatPlot")

data(CeTFdemo)
genes <- unique(c(as.character(NetworkData(CeTFdemo, "network1")[["gene1"]]),
                  as.character(NetworkData(CeTFdemo, "network1")[["gene2"]])))
 

test_that('enrichPlot throws an error when there is none res argument', {
  expect_error(heatPlot(diff = getDE(CeTFdemo, "all"), 
                        showCategory = 10))
})

test_that('enrichPlot throws an error when there is none diff argument', {
  expect_error(heatPlot(res = data.frame(), 
                        showCategory = 10))
})

test_that('enrichPlot throws an error when showCategory argument ir equal 1', {
  expect_error(heatPlot(res = cond1$results, 
                        diff = getDE(CeTFdemo, "all"), 
                        showCategory = 1))
})
