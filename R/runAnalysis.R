#' @title Whole analysis of Regulatory Impact Factors (RIF) and Partial Correlation and Information Theory analysis (PCIT)
#'
#' @description This function uses RIF (Reverter and Chan, 2010) and PCIT (Reverter and Chan, 2008) to run the whole pipeline analysis.
#'
#' @param counts Count data where the rows are genes and coluns the samples (conditions).
#' @param conditions A vector of character identifying the names of conditions (i.e. c("normal", "tumoral"))
#' @param lfc logFoldChange module threshold to define a gene as differentially expressed.
#' @param padj Sifnificance value to define a gene as differentially expressed.
#' @param TFs A vector of character with all transcripts factors of specific organism.
#' @param ncond1 Number of samples that correspond to first condition.
#' @param ncond2 Number of samples that correspond to second condition.
#' @param tolType Type of tolerance given the 3 pairwise correlations (mean, min, max, median).
#' @param diffMethod Method to calculate Differential Expressed (DE) genes (e.g. Reverter and DESEq2)
#'
#' @return A list with two dataframes with the output network of genes/TFs for the first and second conditions, a dataframe with the ket TFs and a dataframe with correspondent separation of genes and TFs. This outputs can be used to generate the networks in Cytoscape.
#'
#' @examples
#' data("simCounts")
#' out <- runAnalysis(counts = simCounts,
#' conditions=c("cond1", "cond2"),
#' lfc = 2.57,
#' padj = 0.05,
#' TFs = paste0("TF_", 1:1000),
#' ncond1 = 10,
#' ncond2= 10,
#' tolType = "mean",
#' diffMethod = "Reverter")
#'
#' @importFrom reshape2 melt
#' @importFrom stats cor
#' @importFrom crayon green
#'
#' @export
runAnalysis <- function(counts, conditions=NULL, lfc = 2.57, padj = 0.05, TFs = NULL, ncond1 = NULL, ncond2= NULL, tolType = "mean", diffMethod = "Reverter") {
  if (!is.data.frame(counts) & !is.matrix(counts)) {
    stop("counts must be a dataframe or a matrix")
  }

  if (length(conditions) != 2) {
    stop("you must input two conditions")
  }

  if (length(TFs) == 0 | !is.character(TFs)) {
    stop("the transcript factors must be a character")
  }

  if (!is.numeric(ncond1) | !is.numeric(ncond2)) {
    stop("the number of conditions must be a numeric greater than zero")
  }


  cat(green(
    "##### STEP 1: TPM filter #####" %+% '\n'
  ))

  colnames(counts)[1:ncond1] <- paste0(colnames(counts)[1:ncond1], "_", conditions[1])
  colnames(counts)[(ncond1+1):(ncond1+ncond2)] <- paste0(colnames(counts)[(ncond1+1):(ncond1+ncond2)], "_", conditions[2])

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
    "##### STEP 2: Differential Expression #####" %+% '\n'
  ))


  if (diffMethod == "Reverter") {
    de.j <- data.frame(cond1 = apply(Clean_Dat[, grep(conditions[1], colnames(Clean_Dat))], 1, mean),
                       cond2 = apply(Clean_Dat[, grep(conditions[2], colnames(Clean_Dat))], 1, mean))
    tmp1 <- de.j[,1] - de.j[,2]
    var=(sum(tmp1 ^ 2)-(sum(tmp1)*sum(tmp1))/(length(tmp1)))/(length(tmp1)-1)
    de.j <- cbind(de.j, diff = (tmp1-(sum(tmp1)/length(tmp1)))/sqrt(var))
    DE_unique <- subset(de.j, abs(de.j$diff) > lfc)
    Target <- rownames(DE_unique)
  } else if (diffMethod == "DESeq2") {
    tmp1 <- counts[rownames(Clean_Dat), ]
    anno <- data.frame(cond = c(rep(conditions[1], ncond1), rep(conditions[2], ncond2)),
                       row.names = colnames(counts))
    Target <- expDiff(counts = tmp1,
                      anno = anno,
                      conditions = conditions,
                      lfc = lfc,
                      padj = padj)
  }

  Background <- rownames(Clean_Dat)

  # Get a list of TFs
  TF_unique <- sort(intersect(TFs, Background))


  cat(green(
    "##### STEP 3: Regulatory Impact Factors analysis #####" %+% '\n'
  ))

  RIF_input <- Clean_Dat[c(Target, TF_unique), c(grep(paste0("_", conditions[1]), colnames(Clean_Dat), fixed = T),
                                                 grep(paste0("_", conditions[2]), colnames(Clean_Dat)))]

  #Run RIF
  RIF_out <- RIF(input = RIF_input,
                 nta = nrow(DE_unique),
                 ntf = length(TF_unique),
                 ncond1 = ncond1,
                 ncond2 = ncond2)

  KeyTF <- subset(RIF_out, sqrt(RIF_out$RIF1 ** 2) > 1.96 | sqrt(RIF_out$RIF2 ** 2) > 1.96)


  cat(green(
    "##### STEP 4: Partial Correlation and Information Theory analysis #####" %+% '\n'
  ))

  net.j <- sort(unique(c(as.character(KeyTF$TF), Target)))

  PCIT_input_cond1 <- Clean_Dat[net.j, grep(conditions[1], colnames(Clean_Dat))]
  PCIT_input_cond2 <- Clean_Dat[net.j, grep(conditions[2], colnames(Clean_Dat))]

  # RUN PCIT ...twice!
  PCIT_out_cond1 <- PCIT(PCIT_input_cond1, tolType = tolType)
  PCIT_out_cond2 <- PCIT(PCIT_input_cond2, tolType = tolType)

  # Collect Lineage-specific connections
  PCIT_out <- cbind(PCIT_out_cond1[[1]], PCIT_out_cond2[[1]])
  Network_cond1 <- subset(PCIT_out, PCIT_out[,4] != 0 & PCIT_out[,8] == 0)[, c(1,2)]
  Network_cond2 <- subset(PCIT_out, PCIT_out[,4] == 0 & PCIT_out[,8] != 0)[, c(1,2)]

  # Count connections for each gene in cond1 and in cond2, focussing on Key TFs
  id.j <- c(as.character(subset(PCIT_out_cond1[[1]], PCIT_out_cond1[[1]]$corr2 != 0)[,1]), as.character(subset(PCIT_out_cond1[[1]], PCIT_out_cond1[[1]]$corr2 != 0)[,2]))
  cond1.j <- as.data.frame(table(id.j))
  id.j <- c(as.character(subset(PCIT_out_cond2[[1]], PCIT_out_cond2[[1]]$corr2 != 0)[,1]), as.character(subset(PCIT_out_cond2[[1]], PCIT_out_cond2[[1]]$corr2 != 0)[,2]))
  cond2.j <- as.data.frame(table(id.j))

  KeyTF_Conn_cond1_cond2 <- merge(KeyTF, merge(cond1.j, cond2.j, by = "id.j"), by.x = "TF", by.y = "id.j")
  KeyTF_Conn_cond1_cond2 <- cbind(KeyTF_Conn_cond1_cond2, KeyTF_Conn_cond1_cond2$Freq.x - KeyTF_Conn_cond1_cond2$Freq.y)
  colnames(KeyTF_Conn_cond1_cond2)[5] <- paste0("freq.", conditions[1])
  colnames(KeyTF_Conn_cond1_cond2)[6] <- paste0("freq.", conditions[2])
  colnames(KeyTF_Conn_cond1_cond2)[7] <- "freq.diff"

  genes <- unique(c(as.character(Network_cond2$gene1), as.character(Network_cond2$gene2), as.character(Network_cond1$gene1), as.character(Network_cond1$gene2)))
  anno <- data.frame(genes = genes,
                     class = ifelse(genes %in% KeyTF$TF, "TF", "gene"))

  return(list(Network_cond1 = Network_cond1,
              Network_cond2 = Network_cond2,
              KeyTF_Conn = KeyTF_Conn_cond1_cond2,
              anno = anno))
}
