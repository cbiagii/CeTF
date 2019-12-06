#' @title ClustCoef
#'
#' @description teste.
#'
#' @param mat teste
#'
#' @return teste.
#'
#' @importFrom stats median
#'
#' @examples
#' teste
#'
#' @export
clustCoef <- function(mat) {
  if (!is.data.frame(mat) & !is.matrix(mat)) {
    stop("mat must be a dataframe or a matrix")
  }

  nGenes <- as.integer(nrow(mat))

  cc <- E <- k <- rep(0, length = nrow(mat))
  idx <- 1:nrow(mat)

  for (i in 1:nGenes) {
    neighbours <- (mat[i,1:nGenes] != 0 | mat[1:nGenes,i] != 0)
    k[i] = sum(neighbours)
    if (k[i] < 2) {
      next()
    }
    E[i] = sum(mat[idx[which(neighbours == T)], idx[which(neighbours == T)]])
  }
  cc = E / ( k * (k-1) )
  return(cc)
}
