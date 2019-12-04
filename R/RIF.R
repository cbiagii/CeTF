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
#' @importFrom crayon green
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

  pb <- txtProgressBar(min = 1, max = ntf, style = 3)
  tmp1 <- NULL
  for (i in 1:ntf) {
    rif1 <- 0
    rif2 <- 0

    for(j in 1:nta) {
      gene_ccorr <- cor(tf[i,1:ncond1], ta[j, 1:ncond1]) #cond1
      if (is.na(gene_ccorr)) { gene_ccorr <- 0 }

      gene_ncorr <- cor(tf[i,(ncond1+1):(ncond1+ncond2)], ta[j, (ncond1+1):(ncond1+ncond2)]) #cond2
      if (is.na(gene_ncorr)) { gene_ncorr <- 0 }

      ave <- (sum(ta[j,1:ncond1])/ncond1 + sum(ta[j,(ncond1+1):(ncond1+ncond2)])/ncond2)/2
      de <- sum(ta[j, 1:ncond1])/ncond1 - sum(ta[j,(ncond1+1):(ncond1+ncond2)])/ncond2
      dw <- gene_ccorr - gene_ncorr

      rif1 = rif1 + ave * de * (dw^2)

      er1 <- sum(ta[j,1:ncond1]/ncond1 * gene_ccorr)
      er2 <- sum(ta[j,(ncond1+1):(ncond1+ncond2)]/ncond2 * gene_ncorr)

      rif2 = rif2 + er1^2 - er2^2
    }
    rif1 = rif1 / nta
    rif2 = rif2 / nta

    tmp1 <- rbind(tmp1, c(rownames(tf)[i], rif1, rif2))
    setTxtProgressBar(pb, i)
  }
  close(pb)

  rownames(tmp1) <- tmp1[,1]
  tmp1 <- tmp1[,-c(1)]
  mode(tmp1) = "numeric"
  tmp1 <- as.data.frame(tmp1)

  rif1 <- (tmp1$V1-mean(tmp1$V1))/sd(tmp1$V1)
  rif2 <- (tmp1$V2-mean(tmp1$V2))/sd(tmp1$V2)

  out <- data.frame(gene = rownames(tmp1),
                    RIF1 = rif1,
                    RIF2 = rif2)

  return(out)
}
