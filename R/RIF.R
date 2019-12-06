#' @title Regulatory Impactor Factors analysis
#'
#' @description teste.
#'
#' @param input teste.
#' @param nta teste.
#' @param ntf teste.
#' @param ncond1 teste.
#' @param ncond2 teste.
#'
#' @return teste.
#'
#' @importFrom utils txtProgressBar
#' @importFrom stats cor
#' @importFrom crayon green %+%
#' @import pbapply pbapply
#'
#' @examples
#' teste
#'
#' @export
RIF <- function(input, nta, ntf, ncond1, ncond2) {
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
