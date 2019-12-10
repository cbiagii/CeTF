#' @title Calculate the clustering coefficient as a percentage
#'
#' @description Given an adjacency matrix, calculate the clustering coefficient as a percentage of non-zero adjacencies.
#'
#' @param mat An adjacency matrix. Calculating the clustering coefficient percentage only makes sense if some connections are zero i.e. no connection.
#'
#' @return A numerical between 0 and 100.
#'
#' @importFrom stats median
#'
#' @examples
#' data(simNorm)
#' clustCoefPercentage(simNorm)
#'
#' @export
clustCoefPercentage <- function(mat) {
    if (!is.data.frame(mat) & !is.matrix(mat)) {
        stop("mat must be a dataframe or a matrix")
    }
    
    tmp <- mat[upper.tri(mat)]
    tmp1 <- tmp[which(tmp != 0)]
    return(length(tmp1)/length(tmp) * 100)
}
