#' @title Smear plot for Differentially Expressed genes and TFs
#'
#' @description Generate an plot for Differentially Expressed (DE) genes and for 
#' specific TF that shows the relationship between log(baseMean) and Difference 
#' of Expression or log2FoldChange. This plot enables to visualize the 
#' distribution of DE genes and TF in both conditions.
#' @param object CeTF class object resulted from \code{\link{runAnalysis}} 
#' function.
#' @param diffMethod Method used to calculate Differentially Expressed (DE)
#' genes (see \code{\link{expDiff}}).
#' @param lfc logFoldChange module threshold to define a gene as differentially 
#' expressed (default: 1.5).
#' @param conditions A vector of characters identifying the names of conditions 
#' (i.e. c('normal', 'tumoral')).
#' @param TF Specify a single TF to be visualized for (used only if type argument 
#' equals TF).
#' @param padjust Significance value to define a gene as differentially expressed 
#' in DESeq2 diffMethod option (default: 0.05).
#' @param label If label is TRUE, shows the names of single TF and its 
#' respectives (default: FALSE).
#' @param type Type of plot (DE or TF). If DE, will plot the smear plot for all 
#' differentally expressed genes and TFs for both conditions, and if TF, will 
#' plot the smear plot for a specific TF and their targets.
#' targets for both conditions (default: FALSE).
#'
#' @return Returns an plot of log2(baseMean) by log2FoldChange or difference of
#' expression for genes and TFs differentially expressed or for a single TF and
#' its targets for both conditions.
#'
#' @importFrom tidyr %>%
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot geom_point theme_bw ylab ggtitle theme geom_hline scale_colour_manual labs
#' @importFrom ggrepel geom_label_repel
#' @importFrom stats na.omit
#'
#'
#' @examples
#' # load the CeTF class object resulted from runAnalysis function
#' data(CeTFdemo)
#'
#' #performing SmearPlot for DE genes and TFs
#' SmearPlot(object = CeTFdemo, 
#'           diffMethod = 'Reverter', 
#'           lfc = 1.5, 
#'           conditions = c('untrt', 'trt'), 
#'           type = 'DE')
#'
#' #performing SmearPlot for DE genes and TFs
#' SmearPlot(object = CeTFdemo,
#'           diffMethod = 'Reverter',
#'           lfc = 1.5,
#'           conditions = c('untrt', 'trt'),
#'           TF = 'ENSG00000205189',
#'           label = FALSE, 
#'           type = 'TF')
#'
#' @export
SmearPlot <- function(object, diffMethod, lfc = 1.5, conditions, TF = NULL, 
    padjust = 0.05, label = FALSE, type = NULL) {
    if (!is(object, "CeTF")) {
        stop("the input must be a CeTF class object")
    }
    if (missing(object)) {
        stop("No \"object\" parameter provided")
    }
    if (missing(diffMethod)) {
        stop("No \"diffMethod\" parameter provided")
    }
    if (missing(conditions)) {
        stop("you must input two conditions")
    }
    if (type == "DE" & !is.null(TF)) {
        stop("the TF is not necessary to DE option")
    }
    
    var1 <- switch(diffMethod, Reverter = "diff", DESeq2 = "log2FoldChange")
    var2 <- switch(diffMethod, Reverter = "Difference of Expression", DESeq2 = "log2FoldChange")
    var3 <- switch(diffMethod, Reverter = "Average Expression", DESeq2 = "log2(baseMean)")
    
    if (type == "DE") {
        tab <- getDE(object, "all")
        tab[["genes"]] <- rownames(tab)
        tab[["Type"]] <- "Not Different"
        
        if (diffMethod == "Reverter") {
            tab[which(tab[[var1]] > lfc), "Type"] <- conditions[1]
            tab[which(tab[[var1]] < -lfc), "Type"] <- conditions[2]
            
            tab[["Type"]][which(tab[["genes"]] %in% NetworkData(object, 
                "tfs") & tab[[var1]] > lfc)] <- paste("TF", conditions[1])
            tab[["Type"]][which(tab[["genes"]] %in% NetworkData(object, 
                "tfs") & tab[[var1]] < -lfc)] <- paste("TF", conditions[2])
            tab[["Type"]][which(tab[["genes"]] %in% NetworkData(object, 
                "tfs") & (tab[[var1]] > -lfc & tab[[var1]] < lfc))] <- "Not Different TF"
            tab[["baseMean"]] <- rowMeans(tab[, c(1, 2)])
        } else if (diffMethod == "DESeq2") {
            tab[which(tab[[var1]] > lfc & tab[["padj"]] < padjust), "Type"] <- conditions[1]
            tab[which(tab[[var1]] < -lfc & tab[["padj"]] < padjust), "Type"] <- conditions[2]
            tab[["Type"]][which(tab[["genes"]] %in% NetworkData(object, 
                "tfs") & tab[[var1]] > lfc & tab[["padj"]] < padjust)] <- paste("TF", 
                conditions[1])
            tab[["Type"]][which(tab[["genes"]] %in% NetworkData(object, 
                "tfs") & tab[[var1]] < -lfc & tab[["padj"]] < padjust)] <- paste("TF", 
                conditions[2])
            tab[["Type"]][which(tab[["genes"]] %in% NetworkData(object, 
                "tfs") & (tab[[var1]] > -lfc & tab[[var1]] < lfc) & tab[["padj"]] < 
                padjust)] <- "Not Different TF"
        }
        
        tab[["Type"]] <- factor(tab[["Type"]], levels = c(conditions[1], 
            conditions[2], "Not Different", paste("TF", conditions[2]), 
            paste("TF", conditions[1]), "Not Different TF"))
        tab <- na.omit(tab)
        
        pt <- ggplot(tab) + geom_point(aes(x = log2(.data[["baseMean"]]), 
            y = .data[[var1]], color = .data[["Type"]]))
        
        if (nrow(subset(tab, tab[["Type"]] == conditions[1])) != 0) {
            pt <- pt + geom_point(data = subset(tab, tab[["Type"]] == conditions[1]), 
                aes(x = log2(subset(tab, tab[["Type"]] == conditions[1])[["baseMean"]]), 
                  y = subset(tab, tab[["Type"]] == conditions[1])[[var1]], 
                  colour = conditions[1]), size = 2)
        }
        
        if (nrow(subset(tab, tab[["Type"]] == conditions[2])) != 0) {
            pt <- pt + geom_point(data = subset(tab, tab[["Type"]] == conditions[2]), 
                aes(x = log2(subset(tab, tab[["Type"]] == conditions[2])[["baseMean"]]), 
                  y = subset(tab, tab[["Type"]] == conditions[2])[[var1]], 
                  colour = conditions[2]), size = 2)
        }
        
        if (nrow(subset(tab, tab[["Type"]] == paste("TF", conditions[1]))) != 
            0) {
            pt <- pt + geom_point(data = subset(tab, tab[["Type"]] == paste("TF", 
                conditions[1])), aes(x = log2(subset(tab, tab[["Type"]] == 
                paste("TF", conditions[1]))[["baseMean"]]), y = subset(tab, 
                tab[["Type"]] == paste("TF", conditions[1]))[[var1]], colour = paste("TF", 
                conditions[1])), size = 2)
        }
        if (nrow(subset(tab, tab[["Type"]] == paste("TF", conditions[2]))) != 
            0) {
            pt <- pt + geom_point(data = subset(tab, tab[["Type"]] == paste("TF", 
                conditions[2])), aes(x = log2(subset(tab, tab[["Type"]] == 
                paste("TF", conditions[2]))[["baseMean"]]), y = subset(tab, 
                tab[["Type"]] == paste("TF", conditions[2]))[[var1]], colour = paste("TF", 
                conditions[2])), size = 2)
        }
        if (nrow(subset(tab, tab[["Type"]] == "Not Different TF")) != 0) {
            pt <- pt + geom_point(data = subset(tab, tab[["Type"]] == "Not Different TF"), 
                aes(x = log2(subset(tab, tab[["Type"]] == "Not Different TF")[["baseMean"]]), 
                  y = subset(tab, tab[["Type"]] == "Not Different TF")[[var1]], 
                  colour = "Not Different TF"), size = 2)
        }
        
        pt <- pt + scale_colour_manual(values = c("red2", "#1515ff", "black", 
            "violetred1", "springgreen", "gray50")) + theme_bw() + ylab(var2) + 
            xlab(var3) + theme(legend.position = "top") + geom_hline(yintercept = lfc, 
            linetype = "dashed", color = "black") + geom_hline(yintercept = -lfc, 
            linetype = "dashed", color = "black") + labs(colour = "")
        return(pt)
    } else if (type == "TF") {
        c1 <- c(as.character(NetworkData(object, "network1")[NetworkData(object, 
            "network1")[, "gene1"] %in% TF, "gene2"]), as.character(NetworkData(object, 
            "network1")[NetworkData(object, "network1")[, "gene2"] %in% 
            TF, "gene1"]))
        c2 <- c(as.character(NetworkData(object, "network2")[NetworkData(object, 
            "network2")[, "gene1"] %in% TF, "gene2"]), as.character(NetworkData(object, 
            "network2")[NetworkData(object, "network2")[, "gene2"] %in% 
            TF, "gene1"]))
        
        if ((length(c1) & length(c2)) == 0) {
            stop("No targets were identified for this TF")
        }
        
        tab <- getDE(object, "all")
        tab[["genes"]] <- rownames(tab)
        tab[["Type"]] <- "None"
        
        tab[which(tab[["genes"]] %in% c1), "Type"] <- paste("Target", conditions[1])
        tab[which(tab[["genes"]] %in% c2), "Type"] <- paste("Target", conditions[2])
        tab[which(tab[["genes"]] == TF), "Type"] <- "TF"
        
        if (diffMethod == "Reverter") {
            tab[["baseMean"]] <- rowMeans(tab[, c(1, 2)])
        }
        tab[["Type"]] <- factor(tab[["Type"]], levels = c("TF", paste("Target", 
            conditions[1]), paste("Target", conditions[2]), "None"))
        
        pt <- ggplot(tab) + geom_point(aes(x = log2(.data[["baseMean"]]), 
            y = .data[[var1]], color = .data[["Type"]])) + geom_point(data = subset(tab, 
            tab[["Type"]] == paste("Target", conditions[1])), aes(x = log2(subset(tab, 
            tab[["Type"]] == paste("Target", conditions[1]))[["baseMean"]]), 
            y = subset(tab, tab[["Type"]] == paste("Target", conditions[1]))[[var1]], 
            colour = paste("Target", conditions[1])), size = 2) + geom_point(data = subset(tab, 
            tab[["Type"]] == paste("Target", conditions[2])), aes(x = log2(subset(tab, 
            tab[["Type"]] == paste("Target", conditions[2]))[["baseMean"]]), 
            y = subset(tab, tab[["Type"]] == paste("Target", conditions[2]))[[var1]], 
            colour = paste("Target", conditions[2])), size = 2) + geom_point(data = subset(tab, 
            tab[["Type"]] == "TF"), aes(x = log2(subset(tab, tab[["Type"]] == 
            "TF")[["baseMean"]]), y = subset(tab, tab[["Type"]] == "TF")[[var1]], 
            colour = "TF"), size = 2) + theme_bw() + ylab(var2) + xlab(var3) + 
            ggtitle(paste0("Smear Plot for ", TF, " and its targets")) + 
            theme(legend.position = "top") + labs(colour = "") + scale_colour_manual(values = c("red", 
            "goldenrod3", "forestgreen", "gray80")) + geom_hline(yintercept = lfc, 
            linetype = "dashed", color = "black") + geom_hline(yintercept = -lfc, 
            linetype = "dashed", color = "black")
        
        if (label) {
            idx <- c(which(tab[["Type"]] == paste("Target", conditions[1])), 
                which(tab[["Type"]] == paste("Target", conditions[2])), 
                which(tab[["Type"]] == "TF"))
            dtsub <- tab[idx, ]
            
            return(pt + geom_label_repel(data = dtsub, aes(x = log2(.data[["baseMean"]]), 
                y = .data[[var1]], color = .data[["Type"]], alpha = 0.3, 
                label = .data[["genes"]]), nudge_y = 1, segment.size = 0.2, 
                segment.color = "grey50", direction = "x", show.legend = FALSE))
        } else {
            return(pt)
        }
    }
}
