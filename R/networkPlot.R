#' @title Final network plot
#'
#' @description Generate the network that integrates the network of a condition with groupGO results.
#'
#' @param netCond Network of a specific condition. Can be found in result of \link[pcitRif]{runAnalysis} (step4 -> network_cond1 or network_cond2)
#' @param netGO Dataframe with the results of \link[pcitRif]{getGroupGO} (second element of list). This result can be decreased applying filters for the pathways selection.
#' @param keyTFs TFs identified as importants by \link[pcitRif]{runAnalysis} (step4 -> keytf)
#'
#' @return A network.
#'
#' @importFrom network network %v% network.vertex.names
#' @importFrom GGally ggnet2
#'
#' @examples
#' \dontrun{
#' #See vignette for more details
#' networkPlot(netCond,
#' netGO,
#' keyTFs)
#' }
#'
#'
#'
#' @export
networkPlot <- function(netCond, netGO, keyTFs) {
  network <- rbind(netCond, netGO)
  net <- network(network, directed = F)

  values <- unique(c(as.character(network$gene1), as.character(network$gene2)))
  pathways <- unique(as.character(netGO$gene1))

  TFs <- as.character(keyTFs$TF)
  gns <- setdiff(values, c(TFs, pathways))

  x <- network.vertex.names(net)
  x <- factor(ifelse(x %in% TFs, "TFs", ifelse(x %in% gns, "Gene", ifelse(x %in% pathways, "GO:BP", ""))))
  net %v% "color" = as.character(x)
  y <- c("#4DAF4A", "#E41A1C", "#377EB8")
  names(y) = levels(x)
  labels = c(pathways, as.character(head(keyTFs, 2)[,1]), as.character(tail(keyTFs, 2)[,1]))

  pt <- ggnet2(net, color="color", color.legend = "", palette=y,
               edge.size = 0.5, edge.color="gray70", label.size=1, alpha = 0.75, size = "degree",
               edge.alpha = 0.5, label = labels, legend.position = "bottom",
               mode = "random") + coord_equal() + guides(size = F)

  return(pt)
}
