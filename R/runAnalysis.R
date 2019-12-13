#' @title Whole analysis of Regulatory Impact Factors (RIF) and Partial Correlation and Information Theory analysis (PCIT)
#'
#' @description This function uses RIF (Reverter and Chan, 2010) and PCIT (Reverter and Chan, 2008) to run the whole pipeline analysis.
#'
#' @param mat Count data where the rows are genes and coluns the samples (conditions).
#' @param conditions A vector of character identifying the names of conditions (i.e. c('normal', 'tumoral')).
#' @param lfc logFoldChange module threshold to define a gene as differentially expressed.
#' @param padj Sifnificance value to define a gene as differentially expressed.
#' @param TFs A vector of character with all transcripts factors of specific organism.
#' @param ncond1 Number of samples that correspond to first condition.
#' @param ncond2 Number of samples that correspond to second condition.
#' @param tolType Type of tolerance given the 3 pairwise correlations (mean, min, max, median).
#' @param diffMethod Method to calculate Differential Expressed (DE) genes (e.g. Reverter or DESEq2).
#' @param data.type Type of input data. If is expression (FPKM, TPM, etc) or counts.
#'
#' @return Returns a pcitRif object with output variables of each step of analysis. The outputs can be used to generate the networks in Cytoscape.
#'
#' @examples
#' data('simCounts')
#' out <- runAnalysis(mat = simCounts,
#' conditions=c('cond1', 'cond2'),
#' lfc = 2.57,
#' padj = 0.05,
#' TFs = paste0('TF_', 1:1000),
#' ncond1 = 10,
#' ncond2= 10,
#' tolType = 'mean',
#' diffMethod = 'Reverter')
#'
#' @importFrom reshape2 melt
#' @importFrom stats cor
#' @importFrom crayon green
#' @importFrom kableExtra kable kable_styling
#'
#' @export
runAnalysis <- function(mat, conditions = NULL, lfc = 2.57, padj = 0.05, TFs = NULL, ncond1 = NULL,
    ncond2 = NULL, tolType = "mean", diffMethod = "Reverter", data.type = "expression") {

    if (data.type == "counts" & !all(mat == floor(mat))) {
        stop("for data.type = counts you must input a table of integers numbers")
    }
    if (!is.data.frame(mat) & !is.matrix(mat)) {
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

    cat(green("##### STEP 1: Data adjustment #####" %+% "\n"))
    colnames(mat)[seq_len(ncond1)] <- paste0(colnames(mat)[seq_len(ncond1)], "_", conditions[1])
    colnames(mat)[(ncond1 + 1):(ncond1 + ncond2)] <- paste0(colnames(mat)[(ncond1 + 1):(ncond1 +
        ncond2)], "_", conditions[2])

    if (data.type == "counts") {
        # Convert counts to TPM
        tpm.j <- countsToTPM(mat)
        tmp1 <- apply(tpm.j != 0, 1, sum)
        tmp2 <- apply(tpm.j, 1, sum)

        # Count the non-zero samples and average expression values for each gene
        ns_ave <- data.frame(sum = apply(tpm.j != 0, 1, sum), mean = as.numeric(ifelse(tmp1 > 0, tmp2/tmp1, 0)))
        ns_ave[is.na(ns_ave)] <- 0
        st <- bivar.awk(ns_ave)

        # Use only genes above half the average for both
        genesok.j <- sort(rownames(subset(ns_ave, ns_ave$sum >= as.numeric(st[[1]])/2 & ns_ave$mean >=
                                              as.numeric(st[[2]])/2)))

        # Normalization
        tmp1 <- tpm.j[sort(genesok.j[genesok.j %in% rownames(tpm.j)]), ]
        Clean_Dat <- normExp(tmp1)
    } else if (data.type == "expression") {
        tpm.j <- mat

        tmp1 <- apply(tpm.j != 0, 1, sum)
        tmp2 <- apply(tpm.j, 1, sum)

        # Count the non-zero samples and average expression values for each gene
        ns_ave <- data.frame(sum = apply(tpm.j != 0, 1, sum), mean = as.numeric(ifelse(tmp1 > 0, tmp2/tmp1, 0)))
        ns_ave[is.na(ns_ave)] <- 0
        st <- bivar.awk(ns_ave)

        # Use only genes above half the average for both
        genesok.j <- sort(rownames(subset(ns_ave, ns_ave$sum >= as.numeric(st[[1]])/2 & ns_ave$mean >= as.numeric(st[[2]])/2)))
        Clean_Dat <- mat[sort(genesok.j[genesok.j %in% rownames(tpm.j)]), ]
    }

    # storing the results of step1 in a list
    list1 <- list(raw = mat, tpm = tpm.j, selected_genes = genesok.j, norm = Clean_Dat)

    cat(green("##### STEP 2: Differential Expression #####" %+% "\n"))
    anno <- data.frame(cond = c(rep(conditions[1], ncond1), rep(conditions[2], ncond2)),
                       row.names = colnames(mat))
    if (diffMethod == "Reverter") {
        Target <- expDiff(exp = Clean_Dat, anno = anno, conditions = conditions, lfc = lfc, padj = padj, diffMethod = diffMethod)
    } else if (diffMethod == "DESeq2") {
        tmp1 <- mat[rownames(Clean_Dat), ]
        Target <- expDiff(exp = tmp1, anno = anno, conditions = conditions, lfc = lfc, padj = padj, diffMethod = diffMethod)
    }

    Background <- rownames(Clean_Dat)

    # Get a list of TFs
    TF_unique <- sort(intersect(TFs, Background))

    # storing the results of step2 in a list
    list2 <- list(de_genes = Target, background_genes = Background, tf = TF_unique)

    cat(green("##### STEP 3: Regulatory Impact Factors analysis #####" %+% "\n"))
    RIF_input <- Clean_Dat[c(Target, TF_unique), c(grep(paste0("_", conditions[1]), colnames(Clean_Dat),
        fixed = TRUE), grep(paste0("_", conditions[2]), colnames(Clean_Dat)))]

    # Run RIF
    RIF_out <- RIF(input = RIF_input, nta = length(Target), ntf = length(TF_unique), ncond1 = ncond1,
        ncond2 = ncond2)

    KeyTF <- subset(RIF_out, sqrt(RIF_out$RIF1^2) > 1.96 | sqrt(RIF_out$RIF2^2) > 1.96)

    # storing the results of step3 in a list
    list3 <- list(input = RIF_input, out = RIF_out, keytf = KeyTF)

    cat(green("##### STEP 4: Partial Correlation and Information Theory analysis #####" %+% "\n"))
    net.j <- sort(unique(c(as.character(KeyTF$TF), Target)))

    PCIT_input_cond1 <- Clean_Dat[net.j, grep(paste0("_", conditions[1]), colnames(Clean_Dat))]
    PCIT_input_cond2 <- Clean_Dat[net.j, grep(paste0("_", conditions[2]), colnames(Clean_Dat))]

    # RUN PCIT ...twice!
    PCIT_out_cond1 <- PCIT(PCIT_input_cond1, tolType = tolType)
    PCIT_out_cond2 <- PCIT(PCIT_input_cond2, tolType = tolType)

    # Collect Lineage-specific connections
    PCIT_out <- cbind(PCIT_out_cond1[[1]], PCIT_out_cond2[[1]])
    Network_cond1 <- subset(PCIT_out, PCIT_out[, 4] != 0 & PCIT_out[, 8] == 0)[, c(1, 2)]
    Network_cond2 <- subset(PCIT_out, PCIT_out[, 4] == 0 & PCIT_out[, 8] != 0)[, c(1, 2)]

    # Count connections for each gene in cond1 and in cond2, focussing on Key TFs
    id.j <- c(as.character(subset(PCIT_out_cond1[[1]], PCIT_out_cond1[[1]]$corr2 != 0)[, 1]),
        as.character(subset(PCIT_out_cond1[[1]], PCIT_out_cond1[[1]]$corr2 != 0)[, 2]))
    cond1.j <- as.data.frame(table(id.j))
    id.j <- c(as.character(subset(PCIT_out_cond2[[1]], PCIT_out_cond2[[1]]$corr2 != 0)[, 1]),
        as.character(subset(PCIT_out_cond2[[1]], PCIT_out_cond2[[1]]$corr2 != 0)[, 2]))
    cond2.j <- as.data.frame(table(id.j))

    KeyTF_Conn_cond1_cond2 <- merge(KeyTF, merge(cond1.j, cond2.j, by = "id.j"), by.x = "TF",
        by.y = "id.j")
    KeyTF_Conn_cond1_cond2 <- cbind(KeyTF_Conn_cond1_cond2, KeyTF_Conn_cond1_cond2$Freq.x - KeyTF_Conn_cond1_cond2$Freq.y)
    colnames(KeyTF_Conn_cond1_cond2)[5] <- paste0("freq.", conditions[1])
    colnames(KeyTF_Conn_cond1_cond2)[6] <- paste0("freq.", conditions[2])
    colnames(KeyTF_Conn_cond1_cond2)[7] <- "freq.diff"
    KeyTF_Conn_cond1_cond2 <- KeyTF_Conn_cond1_cond2[order(KeyTF_Conn_cond1_cond2$freq.diff, decreasing = T), ]

    genes <- unique(c(as.character(Network_cond2$gene1), as.character(Network_cond2$gene2), as.character(Network_cond1$gene1),
        as.character(Network_cond1$gene2)))
    anno <- data.frame(genes = genes, class = ifelse(genes %in% KeyTF$TF, "TF", "gene"))

    # storing the results of step4 in a list
    list4 <- list(genes = net.j, input_cond1 = PCIT_input_cond1, input_cond2 = PCIT_input_cond2,
        out_cond1 = PCIT_out_cond1, out_cond2 = PCIT_out_cond2, network_cond1 = Network_cond1,
        network_cond2 = Network_cond2, keytf = KeyTF_Conn_cond1_cond2, anno = anno)

    return(new("pcitRif", step1 = list1, step2 = list2, step3 = list3, step4 = list4))
}
