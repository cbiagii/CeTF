#' @title Tolerance level between 3 pairwise correlations
#'
#' @description Calculates the local tolerance for every trio of genes.
#'
#' @param a Interactor 1.
#' @param b Interactor 2.
#' @param c Interactor 3.
<<<<<<< HEAD
#' @param tolType Calculation type for tolerance (mean, min, max, median) (default: "mean").
=======
#' @param tolType Calculation type for tolerance (mean, min, max, median).
>>>>>>> 43614a53fc5fd047595c36314fe49c8a0a0915a2
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
    if(missing(a)){stop("No \"a\" parameter provided")}
    if(missing(b)){stop("No \"b\" parameter provided")}
    if(missing(c)){stop("No \"c\" parameter provided")}

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
