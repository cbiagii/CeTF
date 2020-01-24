context("getEnrich")

data(CeTFdemo)
genes <- unique(c(as.character(NetworkData(CeTFdemo, "network1")[["gene1"]]),
                  as.character(NetworkData(CeTFdemo, "network1")[["gene2"]])))


test_that('getEnrich throws an error when there is wrong organism choice', {
  expect_error(getEnrich(organism = "org", fdrThr = 0.05, 
                         minNum = 5, maxNum = 500))
})

test_that('getEnrich throws an error when there is wrong database choice', {
  expect_error(getEnrich(organism = "hsapiens", 
                         database = "db", fdrThr = 0.05, minNum = 5, 
                         maxNum = 500))
})

test_that('getEnrich throws an error when there isnt genes parameter', {
  expect_error(getEnrich(organism = "hsapiens", 
                         database = "geneontology_Biological_Process", 
                         fdrThr = 0.05, minNum = 5, maxNum = 500))
})

test_that('getEnrich throws an error when there isnt refGene parameter', {
  expect_error(getEnrich(organism = "hsapiens", 
                         database = "geneontology_Biological_Process", 
                         genes = genes, fdrThr = 0.05, minNum = 5, maxNum = 500))
})

test_that('getEnrich throws an error when there is wrong GeneType choice', {
  expect_error(getEnrich(organism = "hsapiens", 
                         database = "geneontology_Biological_Process", 
                         genes = genes, refGene = refGenes$Homo_sapiens$ENSEMBL, 
                         GeneType = "nom", fdrThr = 0.05, minNum = 5, maxNum = 500))
})
