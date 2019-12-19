context("netGOplot")

library(pcitRif)
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

test_that('netGOplot throws an error when groupBy choice is TFs and there is none TFs as input', {
  expect_error(netGOplot(netCond = getNet1(pcitrifExample),
                         resultsGO = t1,
                         netGO = t2,
                         anno = getAnno(pcitrifExample),
                         groupBy = 'TFs',
                         label = TRUE))
})

test_that('netGOplot throws an error when groupBy choice is TFs and there is none TFs as input', {
  expect_error(netGOplot(netCond = getNet1(pcitrifExample),
                         resultsGO = t1,
                         netGO = t2,
                         anno = getAnno(pcitrifExample),
                         groupBy = 'genes',
                         label = TRUE))
})
