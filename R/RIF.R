#' @title
#' Regulatory Impact Factors (RIF) analysis
#'
#' @description
#' The RIF algorithm identify critical transcript factors (TF) from
#' gene expression data.
#'
#' @param input A matrix of expression with differentially expressed genes
#' and transcript factors in rows, and the samples in columns.
#' @param nta Number of Differentially Expressed (DE) genes.
#' @param ntf Number of Transcription Factors (TFs).
#' @param ncond1 Number of samples of condition 1.
#' @param ncond2 Number of samples of condition 2.
#'
#' @details
#' The input matrix must have the rows and columns ordered by the following request:
#' \enumerate{
#'    \item \strong{rows}: DE genes followed by TFs;
#'    \item \strong{columns}: samples of condition1 followed by samples of condition2.
#' }
#'
#' @return A dataframe with the regulatory impact factors metric for
#' each transcript factor.
#'
#' @references
#' REVERTER, Antonio et al. Regulatory impact factors: unraveling the
#' transcriptional regulation of complex traits from expression data.
#' Bioinformatics, v. 26, n. 7, p. 896-904, 2010.
#' \url{https://academic.oup.com/bioinformatics/article/26/7/896/212064}
#'
#' @examples
#' # load RIF input example
#' data('RIF_input')
#'
#' # performing RIF analysis
#' RIF_out <- RIF(input = RIF_input,
#'                nta = 104,
#'                ntf = 50,
#'                ncond1 = 10,
#'                ncond2 = 10)
#'
#' @importFrom utils txtProgressBar
#' @importFrom stats cor
#' @importFrom crayon green %+%
#' @import pbapply pbapply
#'
#' @export
RIF <- function(input, nta = NULL, ntf = NULL, ncond1 = NULL,
    ncond2 = NULL) {
<<<<<<< HEAD
    if(!is.data.frame(input) & !is.matrix(input)){stop("input must be a dataframe or a matrix")}
    if(!is.numeric(nta) | !is.numeric(ntf)){stop("nta and ntf variables must be a numeric greater than zero")}
    if(!is.numeric(ncond1) | !is.numeric(ncond2)){stop("the number of conditions must be a numeric greater than zero")}
    if(missing(nta)){stop("No \"nta\" parameter provided")}
    if(missing(ntf)){stop("No \"ntf\" parameter provided")}
    if(missing(ncond1)){stop("No \"ncond1\" parameter provided")}
    if(missing(ncond2)){stop("No \"ncond2\" parameter provided")}
=======
    if (!is.data.frame(input) & !is.matrix(input)) {
        stop("input must be a dataframe or a matrix")
    }
    if (!is.numeric(nta) | !is.numeric(ntf)) {
        stop("nta and ntf variables must be a numeric greater than zero")
    }
    if (!is.numeric(ncond1) | !is.numeric(ncond2)) {
        stop("the number of conditions must be a numeric greater than zero")
    }
>>>>>>> 43614a53fc5fd047595c36314fe49c8a0a0915a2

    cat(green("## Starting Regulatory Impact Factors analysis ##" %+%
        "\n"))

    ta <- input[seq_len(nta), ]
    tf <- input[(nta + 1):nrow(input), ]

    tmp <- pbapply(tf, 1, function(i) {
        rif1 <- 0
        rif2 <- 0
        tmp1 <- apply(ta, 1, function(j) {
            gene_ccorr <- cor(i[seq_len(ncond1)], j[seq_len(ncond1)])  #cond1
            if (is.na(gene_ccorr)) {
                gene_ccorr <- 0
            }
            gene_ncorr <- cor(i[(ncond1 + 1):(ncond1 +
                ncond2)], j[(ncond1 + 1):(ncond1 +
                ncond2)])  #cond2
            if (is.na(gene_ncorr)) {
                gene_ncorr <- 0
            }
            ave <- (sum(j[seq_len(ncond1)])/ncond1 +
                sum(j[(ncond1 + 1):(ncond1 + ncond2)])/ncond2)/2
            de <- sum(j[seq_len(ncond1)])/ncond1 -
                sum(j[(ncond1 + 1):(ncond1 + ncond2)])/ncond2
            dw <- gene_ccorr - gene_ncorr
            rif1 = rif1 + ave * de * (dw^2)
            er1 <- sum(j[seq_len(ncond1)]/ncond1 *
                gene_ccorr)
            er2 <- sum(j[(ncond1 + 1):(ncond1 + ncond2)]/ncond2 *
                gene_ncorr)
            rif2 = rif2 + er1^2 - er2^2
            list((c(rif1 = rif1, rif2 = rif2)))
        })

        rif1 <- sum(vapply(lapply(lapply(tmp1, `[[`,
            1), `[[`, 1), sum, FUN.VALUE = 0))/nta
        rif2 <- sum(vapply(lapply(lapply(tmp1, `[[`,
            1), `[[`, 2), sum, FUN.VALUE = 0))/nta

        list((c(rif1 = rif1, rif2 = rif2)))
    })

    df <- data.frame(matrix(unlist(tmp), nrow = length(tmp),
        byrow = TRUE))

    out <- data.frame(TF = rownames(tf), avgexpr = rowMeans(tf),
        RIF1 = (df$X1 - mean(df[, 1]))/sd(df[, 1]),
        RIF2 = (df$X2 - mean(df[, 2]))/sd(df[, 2]),
        row.names = NULL, stringsAsFactors = FALSE)

    return(out)
}
