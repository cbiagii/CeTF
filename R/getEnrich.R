#' @title Enrichment analysis for genes of network
#'
#' @description Enrichment analysis of a set of genes derived from the network
#' of any condition using WebGestalt interface in R. Given a vector of genes, 
#' this function will return the enrichment related to the selected database.
#'
#' @param organism WebGestaltR supports 12 organisms. Users can use the function 
#' listOrganism() to check available organisms.
#' @param database The functional categories for the enrichment analysis. Users 
#' can use the function listGeneSet() to check the available functional databases 
#' for the selected organism. Multiple databases in a vector are supported too.
#' @param genes Should be an R vector object containing the interesting gene list.
#' @param refGene Should be an R vector object containing the reference gene list.
#' There is a list with reference genes for 5 organisms in this package (see 
#' \code{\link{refGenes}}).
#' @param GeneType The ID type of the genes and refGene (they must be the same type).
#' Users can use the function listIdType() to check the available gene types.
#' @param fdrMethod Has five FDR methods: holm, hochberg, hommel, bonferroni, BH 
#' and BY (default: BH).
#' @param fdrThr The significant threshold for fdrMethod (default: 0.05).
#' @param minNum Will be exclude the categories with the number of annotated 
#' genes less than minNum for enrichment analysis (default: 5).
#' @param maxNum Will be exclude the categories with the number of annotated 
#' genes larger than maxNum for enrichment analysis (default: 500).
#' 
#' @return
#' Returns an list with the results of the enrichment analysis of the genes and 
#' a network with the database ID (column 1) and the corresponding
#' genes (column 2).
#'
#' @importFrom WebGestaltR WebGestaltR listOrganism listGeneSet listIdType
#' @importFrom pbapply pbapply
#'
#' @examples
#' \dontrun{
#' # load the CeTF class object resulted from runAnalysis function
#' data(CeTFdemo)
#'
#' # getting the genes in network of condition 1
#' genes <- unique(c(as.character(NetworkData(CeTFdemo, 'network1')[, 'gene1']),
#'                  as.character(NetworkData(CeTFdemo, 'network1')[, 'gene2'])))
#'
#' # performing getEnrich analysis
#' cond1 <- getEnrich(organism='hsapiens', database='geneontology_Biological_Process', 
#'                    genes=genes, GeneType='ensembl_gene_id', 
#'                    refGene=refGenes$Homo_sapiens$ENSEMBL, 
#'                    fdrMethod = 'BH', fdrThr = 0.05, minNum = 5, maxNum = 500)
#' }
#'
#' @export
getEnrich <- function(organism, database, genes, refGene, GeneType, fdrMethod = "BH", 
    fdrThr = 0.05, minNum = 5, maxNum = 500) {
    if (!organism %in% listOrganism()) {
        stop("Select a valid organism")
    }
    if (!database %in% listGeneSet()[, 1]) {
        stop("Select a valid database to perform the enrichment")
    }
    if (missing(genes)) {
        stop("No genes provided")
    }
    if (missing(refGene)) {
        stop("No refGene provided")
    }
    if (!GeneType %in% listIdType()) {
        stop("Select a valid GeneType for genes")
    }
    if (!fdrMethod %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", 
        "BY")) {
        stop("Select a valid fdrMethod")
    }
    
    res <- WebGestaltR(enrichMethod = "ORA", isOutput = FALSE, organism = organism, 
        enrichDatabase = database, interestGene = genes, interestGeneType = GeneType, 
        referenceGene = refGene, referenceGeneType = GeneType, fdrMethod = fdrMethod, 
        fdrThr = fdrThr, minNum = minNum, maxNum = maxNum)
    
    if (is.null(res)) {
        stop("None pathway enriched: try to use a different set of genes")
    }
    
    colnames(res)[1] <- "ID"
    colnames(res)[11] <- "geneID"
    
    tmp <- pbapply(res, 1, function(x) {
        temp <- NULL
        pathways1 <- NULL
        temp <- strsplit(x[["geneID"]], ";")
        pathways1 <- as.character(x[["ID"]])
        pathways1 <- rep(pathways1, length(temp[[1]]))
        return(data.frame(pathways = pathways1, gc = temp[[1]]))
    })
    tmp <- do.call(rbind, tmp)
    tmp <- data.frame(gene1 = tmp$pathways, gene2 = tmp$gc)
    
    return(list(results = res, netGO = tmp))
}
