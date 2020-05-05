#' @title Network plot of gene-gene/gene-TFs interactions
#'
#' @description Generate the network plot of gene-gene/gene-TFs interactions for both conditions.
#'
#' @param x CeTF object resulted from \code{\link{runAnalysis}} function.
#'
#' @return Returns the network plot for both conditions.
#'
#' @importFrom ggpubr ggarrange
#' @importFrom network network network.vertex.names '%v%' '%v%<-'
#' @importFrom GGally ggnet2
#' @importFrom utils head tail
#' @importFrom ggplot2 coord_equal guides element_rect
#' @importFrom methods is
#'
#' @examples
#' # loading a simulated counts data
#' data('simCounts')
#'
#' # performing runAnalysis function
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
#' # plotting networks conditions
#' netConditionsPlot(out)
#'
#'
#'
#' @export
netConditionsPlot <- function(x) {
    if (!is(x, "CeTF")) {
        stop("the input must be a CeTF class object")
    }
    
    TFs <- NetworkData(x, "keytfs")[["TF"]]
    mainTFs <- NetworkData(x, "keytfs")[order(NetworkData(x, "keytfs")[["freq.diff"]], 
        decreasing = TRUE), ]
    
    # cond1
    nt1 <- network(NetworkData(x, "network1"), directed = FALSE)
    genes1 <- unique(c(as.character(NetworkData(x, "network1")[["gene1"]]), 
        as.character(NetworkData(x, "network1")[["gene2"]])))
    gns1 <- setdiff(genes1, TFs)
    
    map <- data.frame(pattern = c(TFs, gns1), sb = c(rep("TFs", length(TFs)), 
        rep("Gene", length(gns1))), stringsAsFactors = FALSE)
    og <- data.frame(singleValue = network.vertex.names(nt1), stringsAsFactors = FALSE)
    cvt <- vapply(map$pattern, grepl, logical(nrow(og)), og$singleValue)
    og$Category <- map$sb[max.col(cvt, "last")]
    x1 <- factor(og$Category)
    
    nt1 %v% "color" <- as.character(x1)
    
    if (length(levels(x1)) == 1) {
        y = "Set1"
    } else {
        y <- c("#4DAF4A", "#E41A1C")
        names(y) = levels(x1)
    }
    
    labels <- c(as.character(head(mainTFs, 2)[, "TF"]), as.character(tail(mainTFs, 
        2)[, "TF"]))
    pt1 <- ggnet2(nt1, color = "color", color.legend = "", palette = y, 
        edge.size = 0.5, edge.color = "gray70", label.size = 0.5, alpha = 0.75, 
        size = "degree", edge.alpha = 0.5, label = labels, legend.position = "bottom") + 
        coord_equal() + guides(size = FALSE) + ggtitle(paste("Network:", 
        gsub("freq.", "", colnames(NetworkData(x, "keytfs"))[5]))) + theme(panel.background = element_rect(fill = "white", 
        colour = "grey50"))
    
    # cond2
    nt2 <- network(NetworkData(x, "network2"), directed = FALSE)
    genes2 <- unique(c(as.character(NetworkData(x, "network2")[["gene1"]]), 
        as.character(NetworkData(x, "network2")[["gene2"]])))
    gns2 <- setdiff(genes2, TFs)
    
    map <- data.frame(pattern = c(TFs, gns2), sb = c(rep("TFs", length(TFs)), 
        rep("Gene", length(gns2))), stringsAsFactors = FALSE)
    og <- data.frame(singleValue = network.vertex.names(nt2), stringsAsFactors = FALSE)
    cvt <- vapply(map$pattern, grepl, logical(nrow(og)), og$singleValue)
    og$Category <- map$sb[max.col(cvt, "last")]
    x2 <- factor(og$Category)
    
    nt2 %v% "color" <- as.character(x2)
    
    if (length(levels(x1)) == 1) {
        y = "Set1"
    } else {
        y <- c("#4DAF4A", "#E41A1C")
        names(y) = levels(x2)
    }
    
    labels <- c(as.character(head(mainTFs, 2)[, "TF"]), as.character(tail(mainTFs, 
        2)[, "TF"]))
    pt2 <- ggnet2(nt2, color = "color", color.legend = "", palette = y, 
        edge.size = 0.5, edge.color = "gray70", label.size = 0.5, alpha = 0.75, 
        size = "degree", edge.alpha = 0.5, label = labels, legend.position = "bottom") + 
        coord_equal() + guides(size = FALSE) + ggtitle(paste("Network:", 
        gsub("freq.", "", colnames(NetworkData(x, "keytfs"))[6]))) + theme(panel.background = element_rect(fill = "white", 
        colour = "grey50"))
    
    return(ggarrange(pt1, pt2, ncol = 2, common.legend = TRUE, legend = "bottom"))
}
