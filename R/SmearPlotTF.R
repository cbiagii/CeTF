#' @title Smear plot for specifics Transcription Factors
#'
#' @description Generate an plot for a specific TF. The Targets for this TF
#' will be found based on \code{\link{getNet1}} and \code{\link{getNet2}} networks
#' and will be possible to visualize the specific TF with the respective
#' targets taking a look in the relationship between log(baseMean) and Difference
#' of Expression or log2FoldChange.
#'
#' @param object pcitrif class object resulted from \code{\link{runAnalysis}} function.
#' @param diffMethod Method used to calculate Differentially Expressed (DE)
#' genes (see \code{\link{expDiff}}).
#' @param lfc logFoldChange module threshold to define a gene as differentially expressed (default: 1.5).
#' @param conditions A vector of characters identifying the names of conditions (i.e. c('normal', 'tumoral')).
#' @param TF Specify a single TF to be visualized for.
#' @param label If label is TRUE, shows the names of single TF and its respectives
#' targets for both conditions (default: FALSE).
#'
#' @return Returns an plot of log2(baseMean) by log2FoldChange or difference of
#' expression for a single TF and its targets for both conditions.
#'
#' @importFrom tidyr %>%
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot geom_point theme_bw ylab ggtitle theme geom_hline scale_colour_manual labs
#' @importFrom ggrepel geom_label_repel
#'
#' @examples
#' # load the pcitrif class object resulted from runAnalysis function
#' data(pcitrifExample)
#'
# performing SmearPlotTF plot
#' SmearPlotTF(object = pcitrifExample,
#'             diffMethod = 'Reverter',
#'             lfc = 1.5,
#'             conditions = c('untrt', 'trt'),
#'             TF = 'ENSG00000185917',
#'             label = FALSE)
#'
#'
#'
#' @export
SmearPlotTF <- function(object, diffMethod, lfc = 1.5, conditions, TF,
    label = FALSE) {
    if (!is(object, "pcitrif")) {
        stop("the input must be a pcitrif class object")
    }
    if (missing(object)) {
        stop("No \"object\" parameter provided")
    }
    if (missing(diffMethod)) {
        stop("No \"diffMethod\" parameter provided")
    }
    if (length(conditions) != 2) {
        stop("you must input two conditions")
    }
    if (missing(conditions)) {
        stop("No \"conditions\" parameter provided")
    }
    if (missing(TF)) {
        stop("No \"TF\" parameter provided")
    }
    if (length(TF) != 1 | !is.character(TF)) {
        stop("the transcript factor must be a single character")
    }

    var1 <- switch(diffMethod, Reverter = "diff", DESeq2 = "log2FoldChange")
    var2 <- switch(diffMethod, Reverter = "Difference of Expression", DESeq2 = "log2FoldChange")
    var3 <- switch(diffMethod, Reverter = "Average Expression", DESeq2 = "log2(baseMean)")

    c1 <- c(as.character(getNet1(object)[getNet1(object)[, 1] %in% TF,
        2]), as.character(getNet1(object)[getNet1(object)[, 2] %in% TF,
        1]))
    c2 <- c(as.character(getNet2(object)[getNet2(object)[, 1] %in% TF,
        2]), as.character(getNet2(object)[getNet2(object)[, 2] %in% TF,
        1]))

    if ((length(c1) & length(c2)) == 0) {
        stop("No targets were identified for this TF")
    }

    tab <- getDE(object)$DE
    tab$genes <- rownames(tab)
    dt <- tab %>% mutate(Type = ifelse(genes %in% c1, paste("Target", conditions[1]),
        ifelse(genes %in% c2, paste("Target", conditions[2]), "None")))
    dt[which(dt$genes == TF), "Type"] <- "TF"
    if (diffMethod == "Reverter") {
        dt$baseMean <- rowMeans(dt[, c(1, 2)])
    }
    dt$Type <- factor(dt$Type, levels = c("TF", paste("Target", conditions[1]),
        paste("Target", conditions[2]), "None"))

    pt <- ggplot(dt) + geom_point(aes(x = log2(baseMean),
        y = dt[[var1]], color = Type)) + geom_point(data = subset(dt,
        Type == paste("Target", conditions[1])), aes(x = log2(baseMean),
        y = subset(dt, Type == paste("Target", conditions[1]))[[var1]],
        colour = paste("Target", conditions[1])), size = 2) + geom_point(data = subset(dt,
        Type == paste("Target", conditions[2])), aes(x = log2(baseMean),
        y = subset(dt, Type == paste("Target", conditions[2]))[[var1]],
        colour = paste("Target", conditions[2])), size = 2) + geom_point(data = subset(dt,
        Type == "TF"), aes(x = log2(baseMean), y = subset(dt, Type == "TF")[[var1]],
        colour = "TF"), size = 2) + theme_bw() + ylab(var2) + xlab(var3) +
        ggtitle(paste0("Smear Plot for ", TF, " and its targets")) + theme(legend.position = "top") +
        labs(colour = "") + scale_colour_manual(values = c("red", "goldenrod3",
        "forestgreen", "gray80")) + geom_hline(yintercept = lfc, linetype = "dashed",
        color = "black") + geom_hline(yintercept = -lfc, linetype = "dashed",
        color = "black")


    if (label) {
        idx <- c(which(dt$Type == paste("Target", conditions[1])), which(dt$Type ==
            paste("Target", conditions[2])), which(dt$Type == "TF"))
        dtsub <- dt[idx, ]

        return(pt + geom_label_repel(data = dtsub, aes(x = log2(baseMean),
            y = dtsub[[var1]], color = Type, alpha = 0.3, label = genes),
            nudge_y = 1, segment.size = 0.2, segment.color = "grey50",
            direction = "x", show.legend = FALSE))
    } else {
        return(pt)
    }
}
