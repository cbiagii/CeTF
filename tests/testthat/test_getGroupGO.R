context("getGroupGO")

library(org.Hs.eg.db)
data(CeTFdemo)
genes <- unique(c(as.character(NetworkData(CeTFdemo, "network1")[["gene1"]]),
                  as.character(NetworkData(CeTFdemo, "network1")[["gene2"]])))



test_that('getGroupGO throws an error when there are no keyType choice', {
  expect_error(getGroupGO(genes = genes,
                          ont = 'BP',
                          annoPkg = org.Hs.eg.db))
})

test_that('getGroupGO throws an error when there are no annoPkg choice', {
  expect_error(getGroupGO(genes = genes,
                          ont = 'BP',
                          keyType = 'ENSEMBL'))
})
