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
<<<<<<< HEAD
    if(!is.data.frame(mat) & !is.matrix(mat)){stop("mat must be a dataframe or a matrix")}

=======
    if (!is.data.frame(mat) & !is.matrix(mat)) {
        stop("mat must be a dataframe or a matrix")
    }
    
>>>>>>> 43614a53fc5fd047595c36314fe49c8a0a0915a2
    tmp <- mat[upper.tri(mat)]
    tmp1 <- tmp[which(tmp != 0)]
    return(length(tmp1)/length(tmp) * 100)
}
