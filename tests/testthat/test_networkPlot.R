context("networkPlot")

library(org.Hs.eg.db)
data(pcitrifExample)
genes <- unique(c(as.character(getNet1(pcitrifExample)[,1]),
                  as.character(getNet1(pcitrifExample)[,2])))
cond1 <- getGroupGO(genes = genes,
                    ont = 'BP',
                    keyType = 'ENSEMBL',
                    annoPkg = org.Hs.eg.db)
t1 <- head(cond1$results, 12)
t2 <- subset(cond1$netGO, cond1$netGO$gene1 %in% as.character(t1[,1]))

test_that('networkPlot throws an error when there is none netCond input', {
  expect_error(networkPlot(netGO = t2,
                           keyTFs = getKeyTF(pcitrifExample)))
})

test_that('networkPlot throws an error when there is none netGO input', {
  expect_error(networkPlot(netCond = getNet1(pcitrifExample),
                           keyTFs = getKeyTF(pcitrifExample)))
})

test_that('networkPlot throws an error when there is none keyTFs input', {
  expect_error(networkPlot(netCond = getNet1(pcitrifExample),
                           netGO = t2))
})

