#' @title Tolerance
#'
#' @description teste.
#'
#' @param a teste
#' @param b teste
#' @param c teste
#' @param type teste
#'
#' @return teste.
#'
#' @importFrom stats median
#'
#' @examples
#' teste
#'
#' @export
tolerance <- function(a, b, c, type = c("mean", "min", "max", "median")) {
  a_z = (a - b*c) / sqrt((1-b^2) * (1-c^2))
  b_y = (b - a*c) / sqrt((1-a^2) * (1-c^2))
  c_x = (c - a*b) / sqrt((1-a^2) * (1-b^2))

  if (type == "mean") {
    tol <- mean(c(abs(a_z/a), abs(b_y/b), abs(c_x/c)))
  } else if (type == "min") {
    tol <- min(abs(a_z/a), abs(b_y/b), abs(c_x/c))
  } else if (type == "max") {
    tol <- max(abs(a_z/a), abs(b_y/b), abs(c_x/c))
  } else if (type == "median") {
    tol <- median(c(abs(a_z/a), abs(b_y/b), abs(c_x/c)))
  }

  return(tol)
}
