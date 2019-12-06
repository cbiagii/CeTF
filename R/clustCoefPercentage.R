#' @title clustCoefPercentage
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
clustCoefPercentage <- function(mat) {
  if (!is.data.frame(mat) & !is.matrix(mat)) {
    stop("mat must be a dataframe or a matrix")
  }

  tmp <- mat[upper.tri(mat)]
  tmp1 <- tmp[which(tmp != 0)]
  return(length(tmp1)/length(tmp) * 100)
}
