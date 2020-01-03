#' @title Summary statistics from two variables
#'
#' @description Read two columns of data values (say X and Y)
#' and computes summary statistics including N, Mean, SD, Min
#' and Max for X and Y, as well as the correlation between X
#' and Y and the regression of Y on X.
#'
#' @param x A dataframe with two columns (variables).
#'
#' @return Returns an summary statistics for two variables.
#'
#' @importFrom crayon green red blue %+%
#' @importFrom stats lm sd cor
#'
#' @examples
#' # creating a random dataframe with two columns (variables)
#' tab <- data.frame(a = sample(1:1000, 100, replace=TRUE),
#'                   b = sample(1:1000, 100, replace=TRUE))
#'
#' # running bivar.awk function
#' bivar.awk(tab)
#'
#' @export
bivar.awk <- function(x) {
    if (!is.data.frame(x) & !is.matrix(x)) {
        stop("x must be a dataframe or a matrix")
    }
    
    message(green(paste0("## '", colnames(x)[1], "' Column ##") %+% "\n" %+% 
        sprintf("N     =  %s", length(x[, 1])) %+% "\n" %+% sprintf("Mean  =  %s", 
        format(round(mean(x[, 1]), 4), nsmall = 4)) %+% "\n" %+% sprintf("Std   =  %s", 
        format(round(sd(x[, 1]), 4), nsmall = 4)) %+% "\n" %+% sprintf("Min   =  %s", 
        format(round(min(x[, 1]), 4), nsmall = 4)) %+% "\n" %+% sprintf("Max   =  %s", 
        format(round(max(x[, 1]), 4), nsmall = 4)) %+% "\n" %+% "################################" %+% 
        "\n"))
    
    message(red(paste0("## '", colnames(x)[2], "' Column ##") %+% "\n" %+% 
        sprintf("N     =  %s", length(x[, 2])) %+% "\n" %+% sprintf("Mean  =  %s", 
        format(round(mean(x[, 2]), 4), nsmall = 4)) %+% "\n" %+% sprintf("Std   =  %s", 
        format(round(sd(x[, 2]), 4), nsmall = 4)) %+% "\n" %+% sprintf("Min   =  %s", 
        format(round(min(x[, 2]), 4), nsmall = 4)) %+% "\n" %+% sprintf("Max   =  %s", 
        format(round(max(x[, 2]), 4), nsmall = 4)) %+% "\n" %+% "################################" %+% 
        "\n"))
    
    message(blue(sprintf("###" %+% "\n" %+% "Correlation  =  %s", format(round(cor(x[, 
        1], x[, 2]), 7), nsmall = 7)) %+% "\n" %+% sprintf("Regression   =  %s", 
        format(round(lm(x[, 2] ~ x[, 1], data = x)[[1]][[2]], 7), nsmall = 7)) %+% 
        "\n" %+% "###" %+% "\n"))
    
    return(list(mean1 = format(round(mean(x[, 1]), 4), nsmall = 4), mean2 = format(round(mean(x[, 
        2]), 4), nsmall = 4)))
}
