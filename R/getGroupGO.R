#' @title Functional Profile of a gene set at specific GO level
#'
#' @description Functional Profile of a gene set at specific GO level. Given a vector of genes, this function will return the GO profile at a specific level.
#'
#' @param genes Character vector with the genes to perform the functional profile.
#' @param ont One of "MF", "BP", and "CC" subontologies.
#' @param keyType key type of input gene.
#' @param annoPkg Package of annotation of specific organism (OrgDb).
#'
#' @return A list with results of functional profile of genes and a network with the ontology related with which gene.
#'
#' @importFrom clusterProfiler groupGO
#' @importFrom pbapply pbapply
#'
#' @examples
#' \dontrun{
#' results <- getGroupGO(genes,
#' ont = "BP",
#' keyType = "ENSEMBL",
#' annoPkg = org.Hs.eg.db::org.Hs.eg.db)
#' }
#'
#'
#'
#' @export
getGroupGO <- function(genes, ont = "BP", keyType = NULL, annoPkg = NULL) {
  ggo <- groupGO(gene     = as.character(genes),
                 OrgDb    = annoPkg,
                 ont      = ont,
                 readable = F,
                 keyType  = keyType,
                 level    = 3)

  results <- ggo@result
  results <- results[order(results$Count, decreasing = T),]
  results <- results[results$Count>0,]

  tmp <- pbapply(results, 1, function(x){
    temp <- NULL
    pathways1 <- NULL
    temp <- strsplit(x[['geneID']],"/")
    pathways1 <- as.character(x[['ID']])
    pathways1 <- rep(pathways1, length(temp[[1]]))
    return(data.frame(pathways = pathways1,
                      gc = temp[[1]]))
  })
  tmp <- do.call(rbind, tmp)
  tmp <- data.frame(gene1 = tmp$pathways,
                     gene2 = tmp$gc)

  return(list(results = results, netGO = tmp))
}
