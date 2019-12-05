#' @title Partial Correlation and Information Theory analysis
#'
#' @description teste.
#'
#' @param input teste.
#'
#' @return teste.
#'
#' @importFrom reshape2 melt
#' @importFrom stats cor
#' @importFrom crayon green
#'
#' @examples
#' teste
#'
#' @export
PCIT <- function(input){
  "/" <- function(x,y) ifelse(y==0,0,base:::"/"(x,y))

  cat(green(
    "################################" %+% '\n' %+%
      sprintf("Number of genes       =  %s", nrow(input)) %+% '\n' %+%
      sprintf("Number of conditions  =  %s", ncol(input)) %+% '\n' %+%
      "################################" %+% '\n'
  ))

  suppressWarnings(gene_corr <- cor(t(input)))
  gene_corr[is.na(gene_corr)] <- 0

  gene_pcorr <- gene_corr
  gene_pcorr2 <- gene_corr
  for (i in 1:(nrow(gene_pcorr)-2)) {
    if (i%%10==0) { message(paste("Trios for gene", i, sep = "    "))}
    for (j in (i+1):(nrow(gene_pcorr)-1)) {
      for (k in (j+1):nrow(gene_pcorr)) {
        rxy <- gene_pcorr[i,j]
        rxz <- gene_pcorr[i,k]
        ryz <- gene_pcorr[j,k]

        tol <- tolerance(rxy, rxz, ryz, type = "mean")

        if (abs(rxy) < abs(rxz*tol) & abs(rxy) < abs(ryz*tol)) {
          gene_pcorr2[i,j] <- gene_pcorr2[j,i] <- 0
        }

        if (abs(rxz) < abs(rxy*tol) & abs(rxz) < abs(ryz*tol)) {
          gene_pcorr2[i,k] <- gene_pcorr2[k,i] <- 0
        }

        if (abs(ryz) < abs(rxy*tol) & abs(ryz) < abs(rxz*tol)) {
          gene_pcorr2[j,k] <- gene_pcorr2[k,j] <- 0
        }
      }
    }
  }
  rm(gene_pcorr)

  gene_corr <- melt(gene_corr)
  gene_corr <- gene_corr[duplicated(t(apply(gene_corr, 1, sort))), ]
  rownames(gene_corr) <- paste(gene_corr$Var1, gene_corr$Var2, sep = "_")
  gene_corr <- gene_corr[order(gene_corr$Var1), ]

  gene_pcorr2 <- melt(gene_pcorr2)
  rownames(gene_pcorr2) <- paste(gene_pcorr2$Var1, gene_pcorr2$Var2, sep = "_")
  gene_pcorr2 <- gene_pcorr2[rownames(gene_corr), ]

  out <- data.frame(gene1 = gene_corr$Var1,
                    gene2 = gene_corr$Var2,
                    corr1 = round(gene_corr$value, 5),
                    corr2 = round(gene_pcorr2$value, 5))
  rm(gene_corr, gene_corr2)

  return(out)
}
