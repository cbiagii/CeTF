#' @title Regulatory Impact Factors (RIF) analysis
#'
#' @description The RIF algorithm (Reverter and Chan, 2010) identify critical transcript factors (TF) from gene expression data.
#'
#' @param input A matrix of expression with differential expressed genes and transcript factors in rows, and the two conditions in columns.
#' @param nta Number of differential expressed genes.
#' @param ntf Number of transcript factors.
#' @param ncond1 Number of condition 1.
#' @param ncond2 Number of condition 2.
#'
#' @return A dataframe with the regulatory impact factors metric for each transcript factor.
#'
#' @examples
#' data("RIF_input")
#' RIF_out <- RIF(input = RIF_input,
#'nta = 104,
#'ntf = 50,
#'ncond1 = 10,
#'ncond2 = 10)
#'
#' @importFrom utils txtProgressBar
#' @importFrom stats cor
#' @importFrom crayon green %+%
#' @import pbapply pbapply
#'
#' @export
RIF <- function(input, nta=NULL, ntf=NULL, ncond1=NULL, ncond2=NULL) {
  if (!is.data.frame(input) & !is.matrix(input)) {
    stop("input must be a dataframe or a matrix")
  }

  if (!is.numeric(nta) | !is.numeric(ntf)) {
    stop("nta and ntf variables must be a numeric greater than zero")
  }

  if (!is.numeric(ncond1) | !is.numeric(ncond2)) {
    stop("the number of conditions must be a numeric greater than zero")
  }


  cat(green(
    "#################################################" %+% '\n' %+%
      "## Starting Regulatory Impact Factors analysis ##" %+% '\n' %+%
      "#################################################" %+% '\n'
  ))

  ta <- input[1:nta, ]
  tf <- input[(nta+1):nrow(input), ]

  tmp <- pbapply(tf, 1, function(i) {
    rif1 <- 0
    rif2 <- 0
    tmp1 <- apply(ta, 1, function(j) {
      gene_ccorr <- cor(i[1:ncond1], j[1:ncond1]) #cond1
      if (is.na(gene_ccorr)) { gene_ccorr <- 0 }
      gene_ncorr <- cor(i[(ncond1+1):(ncond1+ncond2)], j[(ncond1+1):(ncond1+ncond2)]) #cond2
      if (is.na(gene_ncorr)) { gene_ncorr <- 0 }
      ave <- (sum(j[1:ncond1])/ncond1 + sum(j[(ncond1+1):(ncond1+ncond2)])/ncond2)/2
      de <- sum(j[1:ncond1])/ncond1 - sum(j[(ncond1+1):(ncond1+ncond2)])/ncond2
      dw <- gene_ccorr - gene_ncorr
      rif1 = rif1 + ave * de * (dw^2)
      er1 <- sum(j[1:ncond1]/ncond1 * gene_ccorr)
      er2 <- sum(j[(ncond1+1):(ncond1+ncond2)]/ncond2 * gene_ncorr)
      rif2 = rif2 + er1^2 - er2^2
      list((c(rif1=rif1, rif2=rif2)))
    })

    rif1 <- sum(sapply(lapply(lapply(tmp1, `[[`, 1), `[[`, 1), sum))/nta
    rif2 <- sum(sapply(lapply(lapply(tmp1, `[[`, 1), `[[`, 2), sum))/nta

    list((c(rif1=rif1, rif2=rif2)))
  })

  df <- data.frame(matrix(unlist(tmp), nrow=length(tmp), byrow=T))
  rif1 <- (df$X1 - mean(df[,1]))/sd(df[,1])
  rif2 <- (df$X2 - mean(df[,2]))/sd(df[,2])

  out <- data.frame(gene = rownames(tf),
                     RIF1 = rif1,
                     RIF2 = rif2)

  return(out)
}
