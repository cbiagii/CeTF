context("getGroupGO")

library(org.Hs.eg.db)
data(pcitrifExample)
genes <- unique(c(as.character(getNet1(pcitrifExample)[,1]),
                  as.character(getNet1(pcitrifExample)[,2])))



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
