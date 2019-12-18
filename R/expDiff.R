#' @title Differential expression analysis
#'
#' @description This function returns the differentially expressed genes
#' when comparing two conditions.
#'
#' @param exp Count data where the rows are genes and coluns the samples.
#' @param anno A single column dataframe. The column name must be 'cond',
#' and the rownames must be the names of samples.
#' @param conditions A character vector containing the name of the two
#' conditions. The first name will be selected as reference.
#' @param lfc log2 fold change module threshold to define a gene as
#' differentially expressed (default: 1.5).
#' @param padj Significance value to define a gene as differentially
<<<<<<< HEAD
#' expressed (default: 0.05).
#' @param diffMethod Choose between Reverter or DESeq2 method (default: "Reverter").
=======
#' expressed.
#' @param diffMethod Choose between Reverter or DESeq2 method.
>>>>>>> 43614a53fc5fd047595c36314fe49c8a0a0915a2
#' The DESeq2 method is only for counts data (see details).
#'
#' @return A character with the names of differentially expressed genes.
#'
#' @details
#' The \strong{Reverter} option to diffMethod parameter works as follows:
#' \enumerate{
#'    \item Calculation of mean between samples of each condition for all genes;
#'    \item Subtraction between mean of control condition relative to other condition;
#'    \item Calculation of variance of subtraction previously obtained;
#'    \item The last step calculates the differential expression using the following
#'    formula, where x is the result of substraction (item 2) and var is the
#'    variance calculated in item 3:
#'        \deqn{diff = \frac{x - (sum(x)/length(x))}{\sqrt{var}}}
#' }
#'
#' The \strong{DESeq2} option to diffMethod parameter is recommended only for
#' count data. This method apply the differential expression analysis based
#' on the negative binomial distribution (see \code{\link[DESeq2]{DESeq}}).
#'
#' @examples
#' # loading a simulated counts data
#' data('simCounts')
#'
#' # creating the dataframe with annotation for each sample
#' anno <- data.frame(cond = c(rep('cond1', 10), rep('cond2', 10)),
#'                    row.names = colnames(simCounts))
#'
#' # renaming colums of simulated counts data
#' colnames(simCounts) <- paste(colnames(simCounts), anno$cond, sep = '_')
#'
#' # performing differential expression analysis using Reverter method
#' out <- expDiff(exp = simCounts,
#'                anno = anno,
#'                conditions = c('cond1', 'cond2'),
#'                lfc = 2,
#'                padj = 0.05,
#'                diffMethod = 'Reverter')
#'
#' @references
#' REVERTER, Antonio et al. Simultaneous identification of differential
#' gene expression and connectivity in inflammation, adipogenesis and cancer.
#' Bioinformatics, v. 22, n. 19, p. 2396-2404, 2006.
#' \url{https://academic.oup.com/bioinformatics/article/22/19/2396/240742}
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#' @importFrom stats relevel
#' @importFrom SummarizedExperiment colData
#'
#' @export
expDiff <- function(exp, anno = NULL, conditions = NULL, 
    lfc = 1.5, padj = 0.05, diffMethod = "Reverter") {
<<<<<<< HEAD

    if(missing(anno)){stop("No \"anno\" parameter provided")}
    if(nrow(anno)==0){stop("anno must have some conditions to compare")}
    if(missing(conditions)){stop("No \"conditions\" parameter provided")}
    if(!is.character(conditions) & length(conditions)!=2) {stop("you must input 2 conditions")}
    if(!is.data.frame(exp) & !is.matrix(exp)){stop("exp must be a dataframe or a matrix")}

=======
    if (!is.data.frame(exp) & !is.matrix(exp)) {
        stop("exp must be a dataframe or a matrix")
    }
    
    if (nrow(anno) == 0) {
        stop("anno must have some conditions to compare")
    }
    
    if (!is.character(conditions) & length(conditions) != 
        2) {
        stop("you must input 2 conditions")
    }
    
>>>>>>> 43614a53fc5fd047595c36314fe49c8a0a0915a2
    if (diffMethod == "Reverter") {
        de.j <- data.frame(cond1 = apply(exp[, grep(paste0("_", 
            conditions[1]), colnames(exp))], 1, mean), 
            cond2 = apply(exp[, grep(paste0("_", conditions[2]), 
                colnames(exp))], 1, mean))
        tmp1 <- de.j[, 1] - de.j[, 2]
        var = (sum(tmp1^2) - (sum(tmp1) * sum(tmp1))/(length(tmp1)))/(length(tmp1) - 
            1)
        de.j <- cbind(de.j, diff = (tmp1 - (sum(tmp1)/length(tmp1)))/sqrt(var))
        DE_unique <- rownames(subset(de.j, abs(de.j$diff) > 
            lfc))
    } else if (diffMethod == "DESeq2") {
        ddsHTSeq <- DESeqDataSetFromMatrix(countData = exp, 
            colData = anno, design = ~cond)
        
        colData(ddsHTSeq)[, 1] <- relevel(colData(ddsHTSeq)[, 
            1], ref = conditions[1])
        dds <- DESeq(ddsHTSeq)
        res <- results(dds, contrast = c("cond", conditions[1], 
            conditions[2]))
        DE_unique <- rownames(res)[which(abs(res$log2FoldChange) > 
            lfc & res$padj < padj)]
    }
    return(DE_unique)
}
