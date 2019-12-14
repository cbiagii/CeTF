#' @title Differential expression analysis
#'
#' @description This function calculate the differential expressed genes for two conditions.
#'
#' @param exp Count data where the rows are genes and coluns the samples (conditions).
#' @param anno A dataframe of only one column. The columns must contain the conditions and be named as 'cond', and the rownames must be the names of samples.
#' @param conditions A character vector containing the name of the two conditions. The first name will be selected as reference.
#' @param lfc logFoldChange module threshold to define a gene as differentially expressed.
#' @param padj Sifnificance value to define a gene as differentially expressed.
#' @param diffMethod Choose between Reverter or DESeq2 method. The DESeq2 method is only for counts data.
#'
#' @return A character with the names of differentially expressed genes.
#'
#' @examples
#' data('simCounts')
#' colnames(simCounts) <- paste0('Sample', 1:20)

#' anno <- data.frame(cond = c(rep('cond1', 10), rep('cond2', 10)),
#' row.names = colnames(simCounts))
#'
#' out <- expDiff(exp = simCounts,
#'             anno = anno,
#'             conditions = c('cond1', 'cond2'),
#'             lfc = 2,
#'             padj = 0.05)
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#' @importFrom stats relevel
#' @importFrom SummarizedExperiment colData
#'
#' @export
expDiff <- function(exp, anno = NULL, conditions = NULL, lfc = 1.5, padj = 0.05, diffMethod = "Reverter") {
    if (diffMethod == "Reverter") {
        de.j <- data.frame(cond1 = apply(exp[, grep(paste0("_", conditions[1]), colnames(exp))], 1, mean),
                           cond2 = apply(exp[, grep(paste0("_", conditions[2]), colnames(exp))], 1, mean))
        tmp1 <- de.j[, 1] - de.j[, 2]
        var = (sum(tmp1^2) - (sum(tmp1) * sum(tmp1))/(length(tmp1)))/(length(tmp1) - 1)
        de.j <- cbind(de.j, diff = (tmp1 - (sum(tmp1)/length(tmp1)))/sqrt(var))
        DE_unique <- rownames(subset(de.j, abs(de.j$diff) > lfc))
    } else if (diffMethod == "DESeq2") {
        ddsHTSeq <- DESeqDataSetFromMatrix(countData = exp, colData = anno, design = ~cond)

        colData(ddsHTSeq)[, 1] <- relevel(colData(ddsHTSeq)[, 1], ref = conditions[1])
        dds <- DESeq(ddsHTSeq)
        res <- results(dds, contrast = c("cond", conditions[1], conditions[2]))
        DE_unique <- rownames(res)[which(abs(res$log2FoldChange) > lfc & res$padj < padj)]
    }
    return(DE_unique)
}
