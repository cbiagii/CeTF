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
runAnalysis <- function(counts, conditions=NULL, lfc = 2.57, TFs = NULL, ncond1 = NULL, ncond2= NULL) {

  cat(green(
    "##################################################" %+% '\n' %+%
      "############### STEP 1: TPM filter ###############" %+% '\n' %+%
      "##################################################" %+% '\n'
  ))

  tpm.j <- apply(counts, 2, function(x) {(1000000*x)/sum(x)})
  tmp1 <- apply(tpm.j != 0, 1, sum)
  tmp2 <- apply(tpm.j, 1, sum)

  # Count the non-zero samples and average TPM for each gene
  ns_ave <- data.frame(sum = apply(tpm.j != 0, 1, sum),
                       mean = as.numeric(ifelse(tmp1 > 0, tmp2/tmp1, 0)))
  st <- bivar.awk(ns_ave)

  # Use only genes above half the average for both
  genesok.j <- sort(rownames(subset(ns_ave, ns_ave$sum >= as.numeric(st[[1]])/2 & ns_ave$mean >= as.numeric(st[[2]])/2)))

  # Make meaningful headers (columns 1, 2, 4 and 5 in sample info)
  tmp1 <- tpm.j[sort(genesok.j[genesok.j %in% rownames(tpm.j)]), ]
  Clean_Dat <- apply(tmp1, 2, function(x) {log(x+1)/log(2)})


  cat(green(
    "###################################################" %+% '\n' %+%
      "######### STEP 2: Differential Expression #########" %+% '\n' %+%
      "###################################################" %+% '\n'
  ))

  de.j <- data.frame(cond1 = apply(Clean_Dat[, grep(conditions[1], colnames(Clean_Dat))], 1, mean),
                     cond2 = apply(Clean_Dat[, grep(conditions[2], colnames(Clean_Dat))], 1, mean))

  tmp1 <- de.j[,1] - de.j[,2]

  var=(sum(tmp1 ^ 2)-(sum(tmp1)*sum(tmp1))/(length(tmp1)))/(length(tmp1)-1)
  de.j <- cbind(de.j, diff = (tmp1-(sum(tmp1)/length(tmp1)))/sqrt(var))

  DE_unique <- subset(de.j, abs(de.j$diff) > lfc)
  Background <- rownames(Clean_Dat)
  Target <- rownames(DE_unique)

  # Get a list of TFs
  TF_unique <- sort(intersect(TFs, Background))


  cat(green(
    "##################################################" %+% '\n' %+%
      "### STEP 3: Regulatory Impact Factors analysis ###" %+% '\n' %+%
      "##################################################" %+% '\n'
  ))

  RIF_input <- Clean_Dat[c(Target, TF_unique), c(grep(conditions[1], colnames(Clean_Dat)),
                                                 grep(conditions[2], colnames(Clean_Dat)))]

  #Run RIF
  RIF_out <- RIF(input = RIF_input,
                 nta = nrow(DE_unique),
                 ntf = length(TF_unique),
                 ncond1 = ncond1,
                 ncond2 = ncond2)

  KeyTF <- subset(RIF_out, sqrt(RIF_out$RIF1 ** 2) > 1.96 | sqrt(RIF_out$RIF2 ** 2) > 1.96)


  cat(green(
    "###############################################################" %+% '\n' %+%
      "# STEP 4: Partial Correlation and Information Theory analysis #" %+% '\n' %+%
      "###############################################################" %+% '\n'
  ))

  net.j <- sort(unique(c(as.character(KeyTF$gene), Target)))

  PCIT_input_TT <- Clean_Dat[net.j, grep(conditions[1], colnames(Clean_Dat))]
  PCIT_input_CC <- Clean_Dat[net.j, grep(conditions[2], colnames(Clean_Dat))]

  # RUN PCIT ...twice!
  PCIT_out_TT <- PCIT(PCIT_input_TT)
  PCIT_out_CC <- PCIT(PCIT_input_CC)

  # Collect Lineage-specific connections
  PCIT_out <- cbind(PCIT_out_TT, PCIT_out_CC)
  Network_TT <- subset(PCIT_out, PCIT_out[,4] != 0 & PCIT_out[,8] == 0)[, c(1,2)]
  Network_CC <- subset(PCIT_out, PCIT_out[,4] == 0 & PCIT_out[,8] != 0)[, c(1,2)]

  # Count connections for each gene in TT and in CC, focussing on Key TFs
  id.j <- c(as.character(subset(PCIT_out_TT, PCIT_out_TT$corr2 != 0)[,1]), as.character(subset(PCIT_out_TT, PCIT_out_TT$corr2 != 0)[,2]))
  tt.j <- as.data.frame(table(id.j))
  id.j <- c(as.character(subset(PCIT_out_CC, PCIT_out_CC$corr2 != 0)[,1]), as.character(subset(PCIT_out_CC, PCIT_out_CC$corr2 != 0)[,2]))
  cc.j <- as.data.frame(table(id.j))

  KeyTF_Conn_TT_CC <- merge(KeyTF, merge(tt.j, cc.j, by = "id.j"), by.x = "gene", by.y = "id.j")
  KeyTF_Conn_TT_CC <- KeyTF_Conn_TT_CC[, c(1, 4, 5)]
  KeyTF_Conn_TT_CC <- cbind(KeyTF_Conn_TT_CC, KeyTF_Conn_TT_CC$Freq.x - KeyTF_Conn_TT_CC$Freq.y)
  colnames(KeyTF_Conn_TT_CC) <- c("TF", "freq.Malignant", "freq.Nonmalignant", "freq.diff")

  genes <- unique(c(as.character(Network_CC$gene1), as.character(Network_CC$gene2), as.character(Network_TT$gene1), as.character(Network_TT$gene2)))
  anno <- data.frame(genes = genes,
                     class = ifelse(genes %in% KeyTF$gene, "TF", "gene"))

  return(list(Network_cond1 = Network_TT,
              Network_cond2 = Network_CC,
              KeyTF_Conn = KeyTF_Conn_TT_CC,
              anno = anno))
}
