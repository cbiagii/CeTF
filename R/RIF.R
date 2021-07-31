#' @title
#' Regulatory Impact Factors (RIF) analysis
#'
#' @description
#' The RIF algorithm identify critical transcript factors (TF) from
#' gene expression data.
#'
#' @param input A matrix of expression with differentially expressed genes
#' and transcript factors in rows, and the samples in columns.
#' @param nta Number of Differentially Expressed (DE) genes.
#' @param ntf Number of Transcription Factors (TFs).
#' @param nSamples1 Number of samples of condition 1.
#' @param nSamples2 Number of samples of condition 2.
#'
#' @details
#' The input matrix must have the rows and columns ordered by the following request:
#' \enumerate{
#'    \item \strong{rows}: DE genes followed by TFs;
#'    \item \strong{columns}: samples of condition1 followed by samples of condition2.
#' }
#'
#' @return Returns an dataframe with the regulatory impact factors metric for
#' each transcript factor.
#'
#' @references
#' REVERTER, Antonio et al. Regulatory impact factors: unraveling the
#' transcriptional regulation of complex traits from expression data.
#' Bioinformatics, v. 26, n. 7, p. 896-904, 2010.
#' \url{https://academic.oup.com/bioinformatics/article/26/7/896/212064}
#'
#' @examples
#' # load RIF input example
#' data('RIF_input')
#'
#' # performing RIF analysis
#' RIF_out <- RIF(input = RIF_input,
#'                nta = 104,
#'                ntf = 50,
#'                nSamples1 = 10,
#'                nSamples2 = 10)
#'
#' @importFrom utils txtProgressBar
#' @importFrom stats cor
#'
#' @export
RIF <- function(input, nta = NULL, ntf = NULL, nSamples1 = NULL, nSamples2 = NULL) {
    if (!is.data.frame(input) & !is.matrix(input)) {
        stop("input must be a dataframe or a matrix")
    }
    if (is.null(nta) | is.null(ntf)) {
        stop("nta and ntf variables must be a numeric greater than zero")
    }
    if (is.null(nSamples1) | is.null(nSamples2)) {
        stop("the number of conditions must be a numeric greater than zero")
    }
    
    ta <- input[seq_len(nta), ]
    tf <- input[(nta + 1):nrow(input), ]
    
    mat_ccorr_i <- tf[, seq_len(nSamples1)]
    mat_ccorr_j <- ta[, seq_len(nSamples1)]
    suppressWarnings(mat1 <- cor(t(mat_ccorr_i), t(mat_ccorr_j)))
    mat1[is.na(mat1)] <- 0
    
    mat_ncorr_i <- tf[, (nSamples1 + 1):(nSamples1 + nSamples2)]
    mat_ncorr_j <- ta[, (nSamples1 + 1):(nSamples1 + nSamples2)]
    suppressWarnings(mat2 <- cor(t(mat_ncorr_i), t(mat_ncorr_j)))
    mat2[is.na(mat2)] <- 0
    
    ave <- (rowSums(ta[, seq_len(nSamples1)])/nSamples1 + rowSums(ta[, 
        (nSamples1 + 1):(nSamples1 + nSamples2)])/nSamples2)/2
    
    de <- rowSums(ta[, seq_len(nSamples1)])/nSamples1 - rowSums(ta[, (nSamples1 + 
        1):(nSamples1 + nSamples2)])/nSamples2
    
    tmp <- (mat1 - mat2)
    rif1 <- apply(tmp, 1, function(x) {
        var1 <- sum(ave * de * x^2)
    })
    
    part1 <- rowSums(ta[, seq_len(nSamples1)]/nSamples1)
    er1 <- apply(mat1, 1, function(x) {
        var1 <- part1 * x
    })
    
    part2 <- rowSums(ta[, (nSamples1 + 1):(nSamples1 + nSamples2)]/nSamples2)
    er2 <- apply(mat2, 1, function(x) {
        var1 <- part2 * x
    })
    
    tmp <- er1^2 - er2^2
    
    rif2 <- colSums(tmp)
    
    rif1 <- rif1/nta
    rif2 <- rif2/nta
    
    out <- data.frame(TF = rownames(tf), avgexpr = rowMeans(tf), RIF1 = (rif1 - 
        mean(rif1))/sd(rif1), RIF2 = (rif2 - mean(rif2))/sd(rif2), row.names = NULL, 
        stringsAsFactors = FALSE)
    
    return(out)
}
