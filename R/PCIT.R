#' @title Partial Correlation and Information Theory (PCIT) analysis
#'
#' @description The PCIT algorithm (Reverter and Chan, 2008) identify the significant correlation of a given matrix of expression.
#'
#' @param input A correlation matrix.
#' @param tolType Type of tolerance given the 3 pairwise correlations (mean, min, max, median).
#'
#' @return A list with the significant correlations, raw adjacency matrix and significant adjacency matrix.
#'
#' @examples
#' data('simNorm')
#' results <- PCIT(simNorm)
#' head(results$out, 15)
#'
#' @importFrom reshape2 melt
#' @importFrom stats cor
#' @importFrom crayon green
#'
#' @export
PCIT <- function(input, tolType = "mean") {
    if (!is.data.frame(input) & !is.matrix(input)) {
        stop("input must be a dataframe or a matrix")
    }
    "/" <- function(x, y) ifelse(y == 0, 0, base::"/"(x, y))
    
    cat(green("################################" %+% "\n" %+% sprintf("Number of genes       =  %s", 
        nrow(input)) %+% "\n" %+% sprintf("Number of conditions  =  %s", ncol(input)) %+% 
        "\n" %+% "################################" %+% "\n"))
    
    suppressWarnings(gene_corr <- cor(t(input)))
    gene_corr[is.na(gene_corr)] <- 0
    
    gene_pcorr <- gene_corr
    gene_pcorr2 <- gene_corr
    for (i in seq_len(nrow(gene_pcorr) - 2)) {
        if (i%%10 == 0) {
            message(paste("Trios for gene", i, sep = "    "))
        }
        for (j in (i + 1):(nrow(gene_pcorr) - 1)) {
            for (k in (j + 1):nrow(gene_pcorr)) {
                rxy <- gene_pcorr[i, j]
                rxz <- gene_pcorr[i, k]
                ryz <- gene_pcorr[j, k]
                
                tol <- tolerance(rxy, rxz, ryz, tolType = tolType)
                
                if (abs(rxy) < abs(rxz * tol) & abs(rxy) < abs(ryz * tol)) {
                  gene_pcorr2[i, j] <- gene_pcorr2[j, i] <- 0
                }
                if (abs(rxz) < abs(rxy * tol) & abs(rxz) < abs(ryz * tol)) {
                  gene_pcorr2[i, k] <- gene_pcorr2[k, i] <- 0
                }
                if (abs(ryz) < abs(rxy * tol) & abs(ryz) < abs(rxz * tol)) {
                  gene_pcorr2[j, k] <- gene_pcorr2[k, j] <- 0
                }
            }
        }
    }
    
    gene_corr <- melt(gene_corr)
    gene_corr <- gene_corr[duplicated(t(apply(gene_corr, 1, sort))), ]
    rownames(gene_corr) <- paste(gene_corr$Var1, gene_corr$Var2, sep = "_")
    gene_corr <- gene_corr[order(gene_corr$Var1), ]
    
    tmp1 <- melt(gene_pcorr2)
    rownames(tmp1) <- paste(tmp1$Var1, tmp1$Var2, sep = "_")
    tmp1 <- tmp1[rownames(gene_corr), ]
    
    out <- data.frame(gene1 = gene_corr$Var1, gene2 = gene_corr$Var2, corr1 = round(gene_corr$value, 
        5), corr2 = round(tmp1$value, 5), stringsAsFactors = FALSE)
    
    return(list(tab = out, adj_raw = gene_pcorr, adj_sig = gene_pcorr2))
}
