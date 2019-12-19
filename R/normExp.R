#' @title Normalized expression transformation
#'
#' @description Normalize the expression data of any type of experiment
#' by columns, applying log(x + 1)/log(2).
#'
#' @param tab A matrix or dataframe of expression data (i.e. TPM, counts, FPKM).
#'
#' @return Returns a table with normalized values.
#'
#' @examples
#' # loading a simulated counts data
#' data('simCounts')
#'
#' # getting the TPM matrix from counts
#' tpm <- countsToTPM(simCounts)
#'
#' # normalizing TPM data
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
