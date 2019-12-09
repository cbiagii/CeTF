#' @title Statistics from a data frame
#'
#' @description Provides a report of number of elements, mean, standard deviation, minimum and maximum for each of two columns and calculate the correlation and regression for the columns.
#'
#' @param x Dataframe of two columns.
#'
#' @return Report with statistics.
#'
#' @importFrom crayon green red blue %+%
#' @importFrom stats lm sd cor
#'
#' @examples
#' df <- data.frame(a = sample(1:1000, 100, replace=TRUE), b = sample(1:1000, 100, replace=TRUE))
#' bivar.awk(df)
#'
#' @export
bivar.awk <- function(x) {
  if (!is.data.frame(x) & !is.matrix(x)) {
    stop("x must be a dataframe or a matrix")
  }

  cat(green(
    "########## 1th Column ##########" %+% '\n' %+%
      sprintf("N     =  %s", length(x[,1])) %+% '\n' %+%
      sprintf("Mean  =  %s", format(round(mean(x[,1]), 4), nsmall = 4)) %+% '\n' %+%
      sprintf("Std   =  %s", format(round(sd(x[,1]), 4), nsmall = 4)) %+% '\n' %+%
      sprintf("Min   =  %s", format(round(min(x[,1]), 4), nsmall = 4)) %+% '\n' %+%
      sprintf("Max   =  %s", format(round(max(x[,1]), 4), nsmall = 4)) %+% '\n' %+%
      "################################" %+% '\n'
  ))

  cat(red(
    "########## 2th Column ##########" %+% '\n' %+%
      sprintf("N     =  %s", length(x[,2])) %+% '\n' %+%
      sprintf("Mean  =  %s", format(round(mean(x[,2]), 4), nsmall = 4)) %+% '\n' %+%
      sprintf("Std   =  %s", format(round(sd(x[,2]), 4), nsmall = 4)) %+% '\n' %+%
      sprintf("Min   =  %s", format(round(min(x[,2]), 4), nsmall = 4)) %+% '\n' %+%
      sprintf("Max   =  %s", format(round(max(x[,2]), 4), nsmall = 4)) %+% '\n' %+%
      "################################" %+% '\n'
  ))

  cat(blue(
    sprintf("################################" %+% '\n' %+%
              "Correlation  =  %s", format(round(cor(x[,1], x[,2]), 7), nsmall = 7)) %+% '\n' %+%
      sprintf("Regression   =  %s", format(round(lm(x[,2] ~ x[,1], data = x)[[1]][[2]], 7), nsmall = 7)) %+% '\n' %+%
      "################################" %+% '\n'
  ))

  return(list(mean1 = format(round(mean(x[,1]), 4), nsmall = 4),
              mean2 = format(round(mean(x[,2]), 4), nsmall = 4)))
}
