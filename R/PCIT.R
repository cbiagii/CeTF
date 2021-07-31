#' @title Partial Correlation and Information Theory (PCIT) analysis
#'
#' @description The PCIT algorithm is used for reconstruction of gene co-expression
#' networks (GCN) that combines the concept partial correlation coefficient
#' with information theory to identify significant gene to gene associations
#' defining edges in the reconstruction of GCN.
#'
#' @param input A correlation matrix.
#' @param tolType Type of tolerance (default: 'mean') given the 3 pairwise correlations
#' (see \code{\link{tolerance}}).
#'
#' @return Returns an list with the significant correlations, raw adjacency matrix
#' and significant adjacency matrix.
#'
#' @examples
#' # loading a simulated normalized data
#' data('simNorm')
#'
#' # getting the PCIT results for first 30 genes
#' results <- PCIT(simNorm[1:30, ])
#'
#' # printing PCIT output first 15 rows
#' head(results$tab, 15)
#'
#' @references
#' REVERTER, Antonio; CHAN, Eva KF. Combining partial correlation and
#' an information theory approach to the reversed engineering of gene
#' co-expression networks. Bioinformatics, v. 24, n. 21, p. 2491-2497, 2008.
#' \url{https://academic.oup.com/bioinformatics/article/24/21/2491/192682}
#'
#' @importFrom stats cor
#' @importFrom Matrix Matrix
#'
#' @export
PCIT <- function(input, tolType = "mean") {
    if (!is.data.frame(input) & !is.matrix(input)) {
        stop("input must be a dataframe or a matrix")
    }
    
    "/" <- function(x, y) ifelse(y == 0, 0, base::"/"(x, y))
    
    suppressWarnings(gene_corr <- cor(t(input)))
    gene_corr[is.na(gene_corr)] <- 0
    
    raw_corr <- gene_corr
    
    tt <- switch(tolType, mean = 1, min = 2, max = 3, 1)
    
    mat1 <- Matrix(gene_corr, sparse = TRUE)
    result <- pcitC(cor = mat1, tolType = tt)
    rownames(result) <- colnames(result) <- rownames(mat1)
    
    gene_corr[lower.tri(gene_corr)] = NA
    gene_corr <- data.frame(Var1 = rownames(gene_corr)[row(gene_corr)], 
        Var2 = colnames(gene_corr)[col(gene_corr)], value = c(gene_corr))
    gene_corr <- gene_corr[-which(is.na(gene_corr[["value"]])), ]
    gene_corr <- gene_corr[gene_corr[["value"]] != 1, ]
    rownames(gene_corr) <- paste(gene_corr[["Var1"]], gene_corr[["Var2"]], 
        sep = "_")
    gene_corr <- gene_corr[order(gene_corr[["Var1"]]), ]
    
    tmp1 <- as.matrix(result)
    tmp1[lower.tri(tmp1)] = NA
    tmp1 <- data.frame(Var1 = rownames(tmp1)[row(tmp1)], Var2 = colnames(tmp1)[col(tmp1)], 
        value = c(tmp1))
    tmp1 <- tmp1[-which(is.na(tmp1[["value"]])), ]
    rownames(tmp1) <- paste(tmp1[["Var1"]], tmp1[["Var2"]], sep = "_")
    tmp1 <- tmp1[rownames(gene_corr), ]
    
    out <- data.frame(gene1 = gene_corr[["Var1"]], gene2 = gene_corr[["Var2"]], 
        corr1 = round(gene_corr[["value"]], 5), corr2 = round(tmp1[["value"]], 
            5), stringsAsFactors = FALSE)
    
    return(list(tab = out, adj_raw = raw_corr, adj_sig = as.matrix(result)))
}
