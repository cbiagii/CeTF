#' @title Network plot for each condition
#'
#' @description Generate the network for the conditions.
#'
#' @param x pcitrif object resulted from \link[pcitRif]{runAnalysis} function.
#'
#' @return The network for both conditions.
#'
#' @importFrom ggpubr ggarrange
#' @importFrom network network network.vertex.names %v%
#' @importFrom GGally ggnet2
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
#' singleNetworkPlot(out)
#'
#'
#'
#' @export
singleNetworkPlot <- function(x) {

  TFs <- as.character(x@step4$keytf$TF)
  mainTFs <- x@step4$keytf[order(x@step4$keytf$freq.diff, decreasing = T), ]

  #cond1
  nt1 <- network(x@step4$network_cond1, directed = F)
  genes1 <- unique(c(as.character(x@step4$network_cond1$gene1), as.character(x@step4$network_cond1$gene2)))
  gns1 <- setdiff(genes1, TFs)
  x1 <- network.vertex.names(nt1)
  x1 <- factor(ifelse(x1 %in% TFs, "TFs", ifelse(x1 %in% gns1, "Gene", "")))
  nt1 %v% "color" <- as.character(x1)
  y <- c("#4DAF4A", "#E41A1C")
  names(y) = levels(x1)
  labels <- c(as.character(head(mainTFs, 2)[,1]), as.character(tail(mainTFs, 2)[,1]))
  pt1 <- ggnet2(nt1, color="color", color.legend = "", palette=y, edge.size = 0.5,
                edge.color="gray70", label.size=0.5, alpha = 0.75, size = "degree",
                edge.alpha = 0.5, label = labels, legend.position = "bottom", mode = "random") +
    coord_equal() + guides(size = F) + ggtitle(paste("Network:", gsub("freq.", "", colnames(x@step4$keytf)[5]))) +
    theme(panel.background = element_rect(fill = "white", colour = "grey50"))

  #cond2
  nt2 <- network(x@step4$network_cond2, directed = F)
  genes2 <- unique(c(as.character(x@step4$network_cond2$gene1), as.character(x@step4$network_cond2$gene2)))
  gns2 <- setdiff(genes2, TFs)
  x2 <- network.vertex.names(nt2)
  x2 <- factor(ifelse(x2 %in% TFs, "TFs", ifelse(x2 %in% gns2, "Gene", "")))
  nt2 %v% "color" <- as.character(x2)
  y <- c("#4DAF4A", "#E41A1C")
  names(y) = levels(x2)
  labels <- c(as.character(head(mainTFs, 2)[,1]), as.character(tail(mainTFs, 2)[,1]))
  pt2 <- ggnet2(nt2, color="color", color.legend = "", palette=y, edge.size = 0.5,
                edge.color="gray70", label.size=0.5, alpha = 0.75, size = "degree",
                edge.alpha = 0.5, label = labels, legend.position = "bottom", mode = "random") +
    coord_equal() + guides(size = F) + ggtitle(paste("Network:", gsub("freq.", "", colnames(x@step4$keytf)[6]))) +
    theme(panel.background = element_rect(fill = "white", colour = "grey50"))

  return(ggarrange(pt1, pt2, ncol=2, common.legend = TRUE, legend="bottom"))
}
