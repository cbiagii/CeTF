#' @title Tolerance between 3 pairwise correlations
#'
#' @description Function for calculating the local tolerance given the 3 pairwise correlations between a, b and c.
#'
#' @param a Interactor 1.
#' @param b Interactor 2.
#' @param c Interactor 3.
#' @param tolType Type of tolerance given the 3 pairwise correlations (mean, min, max, median).
#'
#' @return Value of tolerance.
#'
#' @examples
#' tolerance(0.5, -0.65, 0.23, tolType = 'mean')
#' tolerance(0.5, -0.65, 0.23, tolType = 'max')
#' tolerance(0.5, -0.65, 0.23, tolType = 'min')
#' tolerance(0.5, -0.65, 0.23, tolType = 'median')
#'
#' @importFrom stats median
#'
#' @export
tolerance <- function(a, b, c, tolType = "mean") {
    a_z = (a - b * c)/sqrt((1 - b^2) * (1 - c^2))
    b_y = (b - a * c)/sqrt((1 - a^2) * (1 - c^2))
    c_x = (c - a * b)/sqrt((1 - a^2) * (1 - b^2))
    
    if (tolType == "mean") {
        tol <- mean(c(abs(a_z/a), abs(b_y/b), abs(c_x/c)))
    } else if (tolType == "min") {
        tol <- min(abs(a_z/a), abs(b_y/b), abs(c_x/c))
    } else if (tolType == "max") {
        tol <- max(abs(a_z/a), abs(b_y/b), abs(c_x/c))
    } else if (tolType == "median") {
        tol <- median(c(abs(a_z/a), abs(b_y/b), abs(c_x/c)))
    }
    return(tol)
}
