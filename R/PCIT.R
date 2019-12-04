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

        rxy_g_z = (rxy - rxz*ryz) / sqrt((1-rxz^2) * (1-ryz^2))
        rxz_g_y = (rxz - rxy*ryz) / sqrt((1-rxy^2) * (1-ryz^2))
        ryz_g_x = (ryz - rxy*rxz) / sqrt((1-rxy^2) * (1-rxz^2))

        toler = abs(rxy_g_z/rxy)
        toler = toler + abs(rxz_g_y/rxz)
        toler = toler + abs(ryz_g_x/ryz)
        toler = 0.33333 * toler

        if (abs(rxy) < abs(rxz*toler) & abs(rxy) < abs(ryz*toler)) {
          gene_pcorr2[i,j] <- 0
          gene_pcorr2[j,i] <- 0
        }

        if (abs(rxz) < abs(rxy*toler) & abs(rxz) < abs(ryz*toler)) {
          gene_pcorr2[i,k] <- 0
          gene_pcorr2[k,i] <- 0
        }

        if (abs(ryz) < abs(rxy*toler) & abs(ryz) < abs(rxz*toler)) {
          gene_pcorr2[j,k] <- 0
          gene_pcorr2[k,j] <- 0
        }
      }
    }
  }
  rm(gene_pcorr)

  gene_corr <- melt(gene_corr)
  #gene_corr <- subset(gene_corr, gene_corr$value != 1)
  gene_corr <- gene_corr[duplicated(t(apply(gene_corr, 1, sort))), ]
  rownames(gene_corr) <- paste(gene_corr$Var1, gene_corr$Var2, sep = "_")
  gene_corr <- gene_corr[order(gene_corr$Var1), ]

  gene_pcorr2 <- melt(gene_pcorr2)
  #gene_pcorr2 <- subset(gene_pcorr2, gene_pcorr2$value != 1)
  rownames(gene_pcorr2) <- paste(gene_pcorr2$Var1, gene_pcorr2$Var2, sep = "_")
  gene_pcorr2 <- gene_pcorr2[rownames(gene_corr), ]

  out <- data.frame(gene1 = gene_corr$Var1,
                    gene2 = gene_corr$Var2,
                    corr1 = round(gene_corr$value, 5),
                    corr2 = round(gene_pcorr2$value, 5))

  return(out)
}
