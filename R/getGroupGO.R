#' @title Functional Profile of a gene set at specific GO level
#'
#' @description Functional Profile of a gene set at specific GO level.
#' Given a vector of genes, this function will return the GO profile at
#' a specific level.
#'
#' @param genes Character vector with the genes to perform the functional profile.
#' @param ont One of 'MF', 'BP', and 'CC' subontologies (default: 'BP').
#' @param keyType Key type of inputted genes (i.e. 'ENSEMBL', 'SYMBOL', 'ENTREZID').
#' @param annoPkg Package of annotation of specific organism (i.e. org.Hs.eg.db, org.Bt.eg.db, org.Rn.eg.db, etc).
#'
#' @return
#' Returns an list with the results of the functional profile of the genes and a
#' network with the ontologies (column 1) and the corresponding
#' genes (column 2).
#'
#' @importFrom clusterProfiler groupGO
#' @importFrom pbapply pbapply
#'
#' @examples
#' # load the annotation package
#' library(org.Hs.eg.db)
#'
#' # load the pcitrif class object resulted from runAnalysis function
#' data(pcitrifExample)
#'
#' # getting the genes in network of condition 1
#' genes <- unique(c(as.character(getNet1(pcitrifExample)[,1]),
#'                  as.character(getNet1(pcitrifExample)[,2])))
#'
#' # performing getGroupGO analysis
#' cond1 <- getGroupGO(genes = genes,
#'                     ont = 'BP',
#'                     keyType = 'ENSEMBL',
#'                     annoPkg = org.Hs.eg.db)
#'
#'
#'
#' @export
getGroupGO <- function(genes, ont = "BP", keyType = NULL, annoPkg = NULL) {
    if (missing(keyType)) {
        stop("No \"keyType\" parameter provided")
    }
    if (missing(annoPkg)) {
        stop("No \"annoPkg\" parameter provided")
    }
    
    ggo <- groupGO(gene = as.character(genes), OrgDb = annoPkg, ont = ont, readable = FALSE, 
        keyType = keyType, level = 3)
    
    results <- ggo@result
    results <- results[order(results$Count, decreasing = TRUE), ]
    results <- results[results$Count > 0, ]
    
    tmp <- pbapply(results, 1, function(x) {
        temp <- NULL
        pathways1 <- NULL
        temp <- strsplit(x[["geneID"]], "/")
        pathways1 <- as.character(x[["ID"]])
        pathways1 <- rep(pathways1, length(temp[[1]]))
        return(data.frame(pathways = pathways1, gc = temp[[1]]))
    })
    tmp <- do.call(rbind, tmp)
    tmp <- data.frame(gene1 = tmp$pathways, gene2 = tmp$gc)
    
    return(list(results = results, netGO = tmp))
}
