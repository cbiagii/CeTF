#' @title Differential expression analysis
#'
#' @description This function calculate the differential expressed genes for two conditions.
#'
#' @param counts Count data where the rows are genes and coluns the samples (conditions).
#' @param anno A dataframe of only one column. The columns must contain the conditions and be named as 'cond', and the rownames must be the names of samples.
#' @param conditions A character vector containing the name of the two conditions. The first name will be selected as reference.
#' @param lfc logFoldChange module threshold to define a gene as differentially expressed.
#' @param padj Sifnificance value to define a gene as differentially expressed.
#'
#' @return A character with the names of differentially expressed genes.
#'
#' @examples
#' data('simCounts')
#' colnames(simCounts) <- paste0('Sample', 1:20)

#' anno <- data.frame(cond = c(rep('cond1', 10), rep('cond2', 10)),
#' row.names = colnames(simCounts))
#'
#' out <- expDiff(counts = simCounts,
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
expDiff <- function(counts, anno = NULL, conditions = NULL, lfc = 1.5, padj = 0.05) {
    ddsHTSeq <- DESeqDataSetFromMatrix(countData = counts, colData = anno, design = ~cond)
    
    colData(ddsHTSeq)[, 1] <- relevel(colData(ddsHTSeq)[, 1], ref = conditions[1])
    dds <- DESeq(ddsHTSeq)
    res <- results(dds, contrast = c("cond", conditions[1], conditions[2]))
    DE_unique <- rownames(res)[which(abs(res$log2FoldChange) > lfc & res$padj < padj)]
    
    return(DE_unique)
}
