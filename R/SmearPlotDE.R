#' @title Smear plot for Differentially Expressed genes
#'
#' @description Generate an plot for Differentially Expressed (DE) genes that
#' shows the relationship between log(baseMean) and Difference of Expression
#' or log2FoldChange. This plot enables to visualize the distribution of DE genes
#' and TF in both conditions.
#'
#' @param object pcitrif class object resulted from \code{\link{runAnalysis}} function.
#' @param diffMethod Method used to calculate Differentially Expressed (DE)
#' genes: 'Reverter' or 'DESeq2' (see \code{\link{expDiff}}).
#' @param lfc logFoldChange module threshold to define a gene as differentially expressed (default: 1.5).
#' @param padjust Significance value to define a gene as differentially expressed in DESeq2 diffMethod option (default: 0.05).
#' @param conditions A vector of characters identifying the names of conditions (i.e. c('normal', 'tumoral')).
#'
#' @return Returns an plot of log2(baseMean) by log2FoldChange or difference of
#' expression for genes and TF obtained from \code{\link{runAnalysis}} function.
#'
#' @importFrom tidyr %>%
#' @importFrom dplyr mutate
#' @importFrom stats na.omit
#' @importFrom ggplot2 ggplot geom_point theme_bw ylab ggtitle theme geom_hline scale_colour_manual labs
#'
#' @examples
#' # load the pcitrif class object resulted from runAnalysis function
#' data(pcitrifExample)
#'
#' # performing SmearPlotDE plot
#' SmearPlotDE(object = pcitrifExample,
#'             diffMethod = 'Reverter',
#'             lfc = 1.5,
#'             conditions = c('untrt', 'trt'))
#'
#'
#'
#' @export
SmearPlotDE <- function(object, diffMethod, lfc = 1.5, padjust = 0.05,
    conditions) {
    if (!is(object, "pcitrif")) {
        stop("the input must be a pcitrif class object")
    }
    if (missing(object)) {
        stop("No \"object\" parameter provided")
    }
    if (missing(diffMethod)) {
        stop("No \"diffMethod\" parameter provided")
    }
    if (missing(conditions)) {
        stop("No \"conditions\" parameter provided")
    }

    var1 <- switch(diffMethod, Reverter = "diff", DESeq2 = "log2FoldChange")
    var2 <- switch(diffMethod, Reverter = "Difference of Expression", DESeq2 = "log2FoldChange")
    var3 <- switch(diffMethod, Reverter = "Average Expression", DESeq2 = "log2(baseMean)")

    tab <- getDE(object)$DE
    tab$genes <- rownames(tab)

    if (diffMethod == "Reverter") {
        dt <- tab %>% mutate(Type = ifelse(tab[[var1]] > lfc, conditions[1],
            ifelse(tab[[var1]] < -lfc, conditions[2], "Not Different")))
        dt <- dt %>% mutate(Type = ifelse(genes %in% getTF(object) & tab[[var1]] >
            lfc, paste("TF", conditions[1]), ifelse(genes %in% getTF(object) &
            tab[[var1]] < -lfc, paste("TF", conditions[2]), ifelse(genes %in%
            getTF(object) & (tab[[var1]] > -lfc & tab[[var1]] < lfc), "Not Different TF",
            Type))))
        dt$baseMean <- rowMeans(dt[, c(1, 2)])
    } else if (diffMethod == "DESeq2") {
        dt <- tab %>% mutate(Type = ifelse(tab[[var1]] > lfc & padj < padjust,
            conditions[1], ifelse(tab[[var1]] < -lfc & padj < padjust,
                conditions[2], "Not Different")))
        dt <- dt %>% mutate(Type = ifelse((genes %in% getTF(object) & tab[[var1]] >
            lfc & padj < padjust), paste("TF", conditions[1]), ifelse((genes %in%
            getTF(object) & tab[[var1]] < -lfc & padj < padjust), paste("TF",
            conditions[2]), ifelse((genes %in% getTF(object) & (tab[[var1]] >
            -lfc & tab[[var1]] < lfc) & padj < padjust), "Not Different TF",
            Type))))
    }

    dt <- na.omit(dt)

    # creating the different layers for plot
    layer1 <- dt
    layer2 <- subset(dt, Type == conditions[1])
    layer3 <- subset(dt, Type == conditions[2])
    layer4 <- subset(dt, Type == paste("TF", conditions[1]))
    layer5 <- subset(dt, Type == paste("TF", conditions[2]))
    layer6 <- subset(dt, Type == "Not Different TF")

    pt <- ggplot() +
        geom_point(data=layer1, aes(x=log2(baseMean), y=layer1[[var1]]), colour="black", size=2) +
        geom_point(data=layer2, aes(x=log2(baseMean), y=layer2[[var1]]), colour="red2", size=2) +
        geom_point(data=layer3, aes(x=log2(baseMean), y=layer3[[var1]]), colour="#1515ff", size=2) +
        geom_point(data=layer4, aes(x=log2(baseMean), y=layer4[[var1]]), colour="springgreen", size=2) +
        geom_point(data=layer5, aes(x=log2(baseMean), y=layer5[[var1]]), colour="violetred1", size=2) +
        geom_point(data=layer6, aes(x=log2(baseMean), y=layer6[[var1]]), colour="gray50", size=2)+
        theme_bw() + ylab(var2) + xlab(var3) +
        theme(legend.position = "top") +
        geom_hline(yintercept = lfc, linetype = "dashed", color = "black") +
        geom_hline(yintercept = -lfc, linetype = "dashed", color = "black") +
        labs(colour = "")

    return(pt)
}
