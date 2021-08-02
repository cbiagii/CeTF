#' @title Enrichment analysis for genes of network
#'
#' @description Enrichment analysis of a set of genes derived from the network
#' of any condition using clusterProfiler. Given a vector of genes, this function 
#' will return the enrichment related to the selected database.
#'
#' @param genes Should be an R vector object containing the interesting gene list.
#' @param organismDB clusterProfiler supports a lot of different organisms. Users 
#' can check the following link (https://www.bioconductor.org/packages/release/data/annotation/)
#' and search for annotations starting with *org.*.
#' @param keyType The ID type of the input genes (i.e. SYMBOL, ENTREZID, ENSEMBL, etc.).
#' 
#' @param ont The functional categories for the enrichment analysis. The available
#' ontologies are Biological Process (BP), Molecular Function (MF) and 
#' Cellular Component (CC).
#' @param fdrMethod Has five FDR methods: holm, hochberg, hommel, bonferroni, BH,  
#' BY, fdr and none(default: BH).
#' @param fdrThr The significant threshold for selected pathways (default: 0.05).
#' @param minGSSize Will be exclude the categories with the number of annotated 
#' genes less than minGSSize for enrichment analysis (default: 5).
#' @param maxGSSize Will be exclude the categories with the number of annotated 
#' genes larger than maxGSSize for enrichment analysis (default: 500).
#' 
#' @return
#' Returns an list with the results of the enrichment analysis of the genes and 
#' a network with the database ID (column 1) and the corresponding
#' genes (column 2).
#'
#' @importFrom clusterProfiler enrichGO
#' @importFrom dplyr '%>%'
#' 
#' @examples
#' \dontrun{
#' # load the CeTF class object resulted from runAnalysis function
#' library(org.Hs.eg.db)
#' data(CeTFdemo)
#'
#' # getting the genes in network of condition 1
#' genes <- unique(c(as.character(NetworkData(CeTFdemo, 'network1')[, 'gene1']),
#'                  as.character(NetworkData(CeTFdemo, 'network1')[, 'gene2'])))
#'
#' # performing getEnrich analysis
#' cond1 <- getEnrich(genes = genes, organismDB = org.Hs.eg.db, keyType = 'ENSEMBL', 
#'                    ont = 'BP', fdrMethod = "BH", fdrThr = 0.05, minGSSize = 5, 
#'                    maxGSSize = 500)
#' }
#'
#' @export
getEnrich <- function(genes, organismDB, keyType, ont, fdrMethod = "BH", fdrThr = 0.05, 
                      minGSSize = 5, maxGSSize = 500) {
    if (missing(genes)) {
        stop("No genes provided")
    }
    if (!fdrMethod %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", 
                          "BY", "fdr", "none")) {
        stop("Select a valid fdrMethod")
    }
    if (!ont %in% c("BP", "MF", "CC")) {
        stop("Select a valid subontologie: BP, MF or CC")
    }
    
    res <- enrichGO(gene = genes, OrgDb = organismDB, keyType = keyType, 
                    ont = ont, pvalueCutoff = fdrThr, pAdjustMethod = fdrMethod, 
                    minGSSize = minGSSize, maxGSSize = maxGSSize) %>% 
        as.data.frame()
    
    if (is.null(res)) {
        stop("None pathway enriched: try to use a different set of genes")
    }
    
    colnames(res)[1] <- "ID"
    
    tmp <- apply(res, 1, function(x) {
        temp <- NULL
        pathways1 <- NULL
        temp <- strsplit(x[["geneID"]], "\\/")
        pathways1 <- as.character(x[["ID"]])
        pathways1 <- rep(pathways1, length(temp[[1]]))
        return(data.frame(pathways = pathways1, gc = temp[[1]]))
    })
    
    tmp <- do.call(rbind, tmp)
    tmp <- data.frame(gene1 = tmp$pathways, gene2 = tmp$gc)
    
    return(list(results = res, netGO = tmp))
}
