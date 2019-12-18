#' @title Calculate the clustering coefficient as a percentage
#'
#' @description Given an adjacency matrix, calculate the clustering coefficient as a percentage of non-zero adjacencies.
#'
#' @param mat An adjacency matrix. Calculating the clustering coefficient percentage only makes sense if some connections are zero i.e. no connection.
#'
#' @return Returns the clustering coefficient as a porcentage.
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
#' # getting the clustering coefficient as percentage
#' clustCoefPercentage(results$adj_sig)
#'
#' @export
clustCoefPercentage <- function(mat) {
    if(!is.data.frame(mat) & !is.matrix(mat)){stop("mat must be a dataframe or a matrix")}

    tmp <- mat[upper.tri(mat)]
    tmp1 <- tmp[which(tmp != 0)]
    return(length(tmp1)/length(tmp) * 100)
}
