#' @title Whole analysis of Regulatory Impact Factors (RIF) and Partial Correlation and Information Theory analysis (PCIT)
#'
#' @description This function uses RIF and PCIT algorithms to run the whole pipeline analysis.
#' The pipeline is composed by 4 steps:
#' \enumerate{
#'     \item \strong{Step 1:} Data adjustment;
#'     \item \strong{Step 2:} Differential expression analysis;
#'     \item \strong{Step 3:} Regulatory Impact Factors analysis;
#'     \item \strong{Step 4:} Partial Correlation and Information Theory analysis.
#' }
#'
#' @param mat Count data where the rows are genes and coluns the samples (conditions).
#' @param conditions A vector of characters identifying the names of conditions (i.e. c('normal', 'tumor')).
#' @param lfc logFoldChange module threshold to define a gene as differentially expressed (default: 2.57).
#' @param padj Significance value to define a gene as differentially expressed (default: 0.05).
#' @param TFs A vector of character with all transcripts factors of specific organism.
#' @param nSamples1 Number of samples that correspond to first condition.
#' @param nSamples2 Number of samples that correspond to second condition.
#' @param tolType Tolerance calculation type (see \code{\link{tolerance}}) (default: 'mean').
#' @param diffMethod Method to calculate Differentially Expressed (DE) genes (see \code{\link{expDiff}}) (default: 'Reverter')
#' @param data.type Type of input data. If is \emph{expression} (FPKM, TPM, etc) or \emph{counts.}
#'
#' @return Returns an CeTF class object with output variables of each step of analysis.
#'
#' @examples
#' data('simCounts')
#' out <- runAnalysis(mat = simCounts,
#'                    conditions=c('cond1', 'cond2'),
#'                    lfc = 3,
#'                    padj = 0.05,
#'                    TFs = paste0('TF_', 1:1000),
#'                    nSamples1 = 10,
#'                    nSamples2= 10,
#'                    tolType = 'mean',
#'                    diffMethod = 'Reverter',
#'                    data.type = 'counts')
#'
#' @seealso \code{\link{CeTF-class}}
#'
#' @importFrom reshape2 melt
#' @importFrom stats cor
#' @importFrom crayon green
#' @importFrom S4Vectors SimpleList
#'
#' @export
runAnalysis <- function(mat, conditions = NULL, lfc = 2.57, padj = 0.05, 
    TFs = NULL, nSamples1 = NULL, nSamples2 = NULL, tolType = "mean", diffMethod = "Reverter", 
    data.type = NULL) {
    
    if (data.type == "counts" & !all(mat == floor(mat))) {
        stop("for data.type = counts you must input a table of integers numbers")
    }
    if (!is.data.frame(mat) & !is.matrix(mat)) {
        stop("counts must be a dataframe or a matrix")
    }
    if (missing(conditions)) {
        stop("No \"conditions\" parameter provided")
    }
    if (length(conditions) != 2) {
        stop("you must input two conditions")
    }
    if (is.null(TFs)) {
        stop("the transcript factors must be a character")
    }
    if (is.null(nSamples1) | is.null(nSamples2)) {
        stop("the number of conditions must be a numeric greater than zero")
    }
    if (is.null(data.type)) {
        stop("the data type muste be counts or expression")
    }
    
    message(green("##### STEP 1: Data adjustment #####" %+% "\n"))
    mat <- as.matrix(mat)
    colnames(mat)[seq_len(nSamples1)] <- paste0(colnames(mat)[seq_len(nSamples1)], 
        "_", conditions[1])
    colnames(mat)[(nSamples1 + 1):(nSamples1 + nSamples2)] <- paste0(colnames(mat)[(nSamples1 + 
        1):(nSamples1 + nSamples2)], "_", conditions[2])
    
    if (data.type == "counts") {
        # Convert counts to TPM
        tpm.j <- apply(mat, 2, function(x) {
            (1e+06 * x)/sum(x)
        })
        tmp1 <- apply(tpm.j != 0, 1, sum)
        tmp2 <- apply(tpm.j, 1, sum)
        
        # Count the non-zero samples and average expression values for each
        # gene
        ns_ave <- data.frame(sum = apply(tpm.j != 0, 1, sum), mean = as.numeric(ifelse(tmp1 > 
            0, tmp2/tmp1, 0)))
        ns_ave[is.na(ns_ave)] <- 0
        st <- bivar.awk(ns_ave)
        
        # Use only genes above half the average for both
        genesok.j <- sort(rownames(subset(ns_ave, ns_ave$sum >= as.numeric(st[[1]])/2 & 
            ns_ave$mean >= as.numeric(st[[2]])/2)))
        
        # Normalization
        tmp1 <- tpm.j[sort(genesok.j[genesok.j %in% rownames(tpm.j)]), 
            ]
        Clean_Dat <- normExp(tmp1)
    } else if (data.type == "expression") {
        tpm.j <- mat
        
        tmp1 <- apply(tpm.j != 0, 1, sum)
        tmp2 <- apply(tpm.j, 1, sum)
        
        # Count the non-zero samples and average expression values for each
        # gene
        ns_ave <- data.frame(sum = apply(tpm.j != 0, 1, sum), mean = as.numeric(ifelse(tmp1 > 
            0, tmp2/tmp1, 0)))
        ns_ave[is.na(ns_ave)] <- 0
        st <- bivar.awk(ns_ave)
        
        # Use only genes above half the average for both
        genesok.j <- sort(rownames(subset(ns_ave, ns_ave$sum >= as.numeric(st[[1]])/2 & 
            ns_ave$mean >= as.numeric(st[[2]])/2)))
        Clean_Dat <- mat[sort(genesok.j[genesok.j %in% rownames(tpm.j)]), 
            ]
    }
    
    
    message(green("##### STEP 2: Differential Expression #####" %+% "\n"))
    anno <- data.frame(cond = c(rep(conditions[1], nSamples1), rep(conditions[2], 
        nSamples2)), row.names = colnames(mat))
    if (diffMethod == "Reverter") {
        Target <- expDiff(exp = Clean_Dat, anno = anno, conditions = conditions, 
            lfc = lfc, padj = padj, diffMethod = diffMethod)
    } else if (diffMethod == "DESeq2") {
        tmp1 <- mat[rownames(Clean_Dat), ]
        Target <- expDiff(exp = tmp1, anno = anno, conditions = conditions, 
            lfc = lfc, padj = padj, diffMethod = diffMethod)
    }
    
    Background <- rownames(Clean_Dat)
    
    # Get a list of TFs
    TF_unique <- sort(intersect(TFs, Background))
    
    
    message(green("##### STEP 3: Regulatory Impact Factors analysis #####" %+% 
        "\n"))
    RIF_input <- Clean_Dat[c(rownames(Target$DE_unique), TF_unique), c(grep(paste0("_", 
        conditions[1]), colnames(Clean_Dat), fixed = TRUE), grep(paste0("_", 
        conditions[2]), colnames(Clean_Dat)))]
    
    # Run RIF
    RIF_out <- RIF(input = RIF_input, nta = nrow(Target$DE_unique), ntf = length(TF_unique), 
        nSamples1 = nSamples1, nSamples2 = nSamples2)
    
    KeyTF <- subset(RIF_out, sqrt(RIF_out$RIF1^2) > 1.96 | sqrt(RIF_out$RIF2^2) > 
        1.96)
    
    
    message(green("##### STEP 4: Partial Correlation and Information Theory analysis #####" %+% 
        "\n"))
    net.j <- sort(unique(c(as.character(KeyTF$TF), rownames(Target$DE_unique))))
    
    PCIT_input_cond1 <- Clean_Dat[net.j, grep(paste0("_", conditions[1]), 
        colnames(Clean_Dat))]
    PCIT_input_cond2 <- Clean_Dat[net.j, grep(paste0("_", conditions[2]), 
        colnames(Clean_Dat))]
    
    # RUN PCIT ...twice!
    PCIT_out_cond1 <- PCIT(PCIT_input_cond1, tolType = tolType)
    PCIT_out_cond2 <- PCIT(PCIT_input_cond2, tolType = tolType)
    
    # Collect Lineage-specific connections
    PCIT_out <- cbind(PCIT_out_cond1[[1]], PCIT_out_cond2[[1]])
    Network_cond1 <- subset(PCIT_out, PCIT_out[, 4] != 0 & PCIT_out[, 8] == 
        0)[, c(1, 2)]
    Network_cond2 <- subset(PCIT_out, PCIT_out[, 4] == 0 & PCIT_out[, 8] != 
        0)[, c(1, 2)]
    
    # Count connections for each gene in cond1 and in cond2, focussing on
    # Key TFs
    id.j <- c(as.character(subset(PCIT_out_cond1[[1]], PCIT_out_cond1[[1]]$corr2 != 
        0)[, "gene1"]), as.character(subset(PCIT_out_cond1[[1]], PCIT_out_cond1[[1]]$corr2 != 
        0)[, "gene2"]))
    cond1.j <- as.data.frame(table(id.j))
    id.j <- c(as.character(subset(PCIT_out_cond2[[1]], PCIT_out_cond2[[1]]$corr2 != 
        0)[, "gene1"]), as.character(subset(PCIT_out_cond2[[1]], PCIT_out_cond2[[1]]$corr2 != 
        0)[, "gene2"]))
    cond2.j <- as.data.frame(table(id.j))
    
    KeyTF_Conn_cond1_cond2 <- merge(KeyTF, merge(cond1.j, cond2.j, by = "id.j"), 
        by.x = "TF", by.y = "id.j")
    KeyTF_Conn_cond1_cond2 <- cbind(KeyTF_Conn_cond1_cond2, KeyTF_Conn_cond1_cond2$Freq.x - 
        KeyTF_Conn_cond1_cond2$Freq.y)
    colnames(KeyTF_Conn_cond1_cond2)[5] <- paste0("freq.", conditions[1])
    colnames(KeyTF_Conn_cond1_cond2)[6] <- paste0("freq.", conditions[2])
    colnames(KeyTF_Conn_cond1_cond2)[7] <- "freq.diff"
    KeyTF_Conn_cond1_cond2 <- KeyTF_Conn_cond1_cond2[order(KeyTF_Conn_cond1_cond2$freq.diff, 
        decreasing = TRUE), ]
    
    genes <- unique(c(as.character(Network_cond2$gene1), as.character(Network_cond2$gene2), 
        as.character(Network_cond1$gene1), as.character(Network_cond1$gene2)))
    anno <- data.frame(genes = genes, class = ifelse(genes %in% KeyTF$TF, 
        "TF", "gene"))
    
    
    # storing the results in CeTF class object
    output <- new("CeTF")
    output@Data$data <- SimpleList(raw = mat, tpm = tpm.j, norm = Clean_Dat)
    output@DE$data <- SimpleList(DE = Target)
    output@Input$data <- SimpleList(RIF_input = RIF_input, PCIT_input_cond1 = PCIT_input_cond1, 
        PCIT_input_cond2 = PCIT_input_cond2)
    output@Output$data <- SimpleList(RIF_out = RIF_out, PCIT_out_cond1 = PCIT_out_cond1, 
        PCIT_out_cond2 = PCIT_out_cond2)
    output@Network$data <- SimpleList(net_cond1 = Network_cond1, net_cond2 = Network_cond2, 
        keyTF = KeyTF_Conn_cond1_cond2, tfs = TF_unique, anno = anno)
    return(output)
}
