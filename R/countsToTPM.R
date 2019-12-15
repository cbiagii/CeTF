#' @title Counts to TPM
#'
#' @description Convert the count data to TPM
#'
#' @param counts A matrix or dataframe of counts
#'
#' @return Returns a table with TPM values
#'
#' @examples
#' data('simCounts')
#' tpm <- countsToTPM(simCounts)
#'
#'
#'
#' @export
countsToTPM <- function(counts) {
    if (!is.data.frame(counts) & !is.matrix(counts)) {
        stop("input must be a count dataframe or a matrix")
    }
    
    return(apply(counts, 2, function(x) {
        (1e+06 * x)/sum(x)
    }))
}
