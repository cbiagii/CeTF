#' @title Calculate the clustering coefficient
#'
#' @description Calculate the clustering coefficient for an adjacency matrix.
#'
#' @param mat An adjacency matrix. Calculating the clustering coefficient only makes sense if some connections are zero i.e. no connection.
#'
#' @return Returns the clustering coefficient(s) for the adjacency matrix.
#'
#' @importFrom stats median
#'
#' @references
#' Nathan S. Watson-Haigh, Haja N. Kadarmideen, and Antonio Reverter (2010).
#' PCIT: an R package for weighted gene co-expression networks based on
#' partial correlation and information theory approaches. Bioinformatics.
#' 26(3) 411-413. \url{https://academic.oup.com/bioinformatics/article/26/3/411/215002}
#'
#' @examples
#' # loading a simulated counts data
#' data('simNorm')
#'
#' # running PCIT analysis
#' results <- PCIT(simNorm)
#'
#' # getting the clustering coefficient
#' clustCoef(results$adj_sig)
#'
#' @export
clustCoef <- function(mat) {
    if (!is.data.frame(mat) & !is.matrix(mat)) {
        stop("mat must be a dataframe or a matrix")
    }
    
    nGenes <- as.integer(nrow(mat))
    
    cc <- E <- k <- rep(0, length = nrow(mat))
    idx <- seq_len(nrow(mat))
    
    for (i in seq_len(nGenes)) {
        neighbours <- (mat[i, seq_len(nGenes)] != 0 | 
            mat[seq_len(nGenes), i] != 0)
        k[i] = sum(neighbours)
        if (k[i] < 2) {
            (next)()
        }
        E[i] = sum(mat[idx[which(neighbours == TRUE)], 
            idx[which(neighbours == TRUE)]])
    }
    cc = E/(k * (k - 1))
    return(cc)
}
