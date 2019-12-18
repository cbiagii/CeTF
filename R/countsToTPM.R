#' @title Counts data to TPM
#'
#' @description Convert the count data to TPM. Note: this function
#' is written very simply and can be easily altered to produce other
#' behavior by examining the source code.
#'
#' @param counts A matrix or dataframe with counts data, where the rows
#' contains the genes and the columns the samples.
#'
#' @return Returns an table with TPM values.
#'
#' @examples
#' # loading a simulated counts data
#' data('simCounts')
#'
#' # getting the TPM matrix from counts
#' tpm <- countsToTPM(simCounts)
#'
#'
#'
#' @export
countsToTPM <- function(counts) {
    if(!is.data.frame(counts) & !is.matrix(counts)){stop("input must be a count dataframe or a matrix")}

    return(apply(counts, 2, function(x) {
        (1e+06 * x)/sum(x)
    }))
}
