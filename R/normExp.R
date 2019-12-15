#' @title Normalization of data
#'
#' @description Normalize the expression data of any type of experiment
#' by columns, applying log(x + 1)/log(2).
#'
#' @param tab A matrix or dataframe of expression data
#'
#' @return Returns a table with noralized values
#'
#' @examples
#' data('simCounts')
#' tpm <- countsToTPM(simCounts)
#' norm <- normExp(tpm)
#'
#'
#'
#' @export
normExp <- function(tab) {
    if (!is.data.frame(tab) & !is.matrix(tab)) {
        stop("input must be a expression dataframe or a matrix")
    }
    return(apply(tab, 2, function(x) {
        log(x + 1)/log(2)
    }))
}
