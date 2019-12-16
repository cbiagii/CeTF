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
#' @importFrom network network network.vertex.names '%v%' '%v%<-'
#' @importFrom GGally ggnet2
#' @importFrom utils head tail
#' @importFrom ggplot2 coord_equal guides
#'
#' @examples
#' library(org.Hs.eg.db)
#'
#' data("pcitrifExample")
#'
#' genes <- unique(c(as.character(getNet1(pcitrifExample)[,1]),
#'                  as.character(getNet1(pcitrifExample)[,2])))
#'
#' cond1 <- getGroupGO(genes = genes,
#'                     ont = "BP",
#'                     keyType = "ENSEMBL",
#'                     annoPkg = org.Hs.eg.db)
#'
#' t1 <- head(cond1$results, 12)
#' #Subsetting the network for the conditions to make available only the 12 nodes subsetted
#' t2 <- subset(cond1$netGO, cond1$netGO$gene1 %in% as.character(t1[,1]))
#'
#' networkPlot(netCond = getNet1(pcitrifExample),
#'             netGO = t2,
#'             keyTFs = getKeyTF(pcitrifExample))
#'
#'
#'
#' @export
networkPlot <- function(netCond, netGO, keyTFs) {
    network <- rbind(netCond, netGO)
    net <- network(network, directed = FALSE)

    values <- unique(c(as.character(network$gene1),
        as.character(network$gene2)))
    pathways <- unique(as.character(netGO$gene1))

    TFs <- as.character(keyTFs$TF)
    gns <- setdiff(values, c(TFs, pathways))

    x <- network.vertex.names(net)
    x <- factor(ifelse(x %in% TFs, "TFs", ifelse(x %in%
        gns, "Gene", ifelse(x %in% pathways, "GO:BP",
        ""))))
    net %v% "color" = as.character(x)
    y <- c("#4DAF4A", "#E41A1C", "#377EB8")
    names(y) = levels(x)
    labels = c(pathways, as.character(head(keyTFs,
        2)[, 1]), as.character(tail(keyTFs, 2)[, 1]))

    pt <- ggnet2(net, color = "color", color.legend = "",
        palette = y, edge.size = 0.5, edge.color = "gray70",
        label.size = 1, alpha = 0.75, size = "degree",
        edge.alpha = 0.5, label = labels, legend.position = "bottom",
        mode = "random") + coord_equal() + guides(size = FALSE)

    return(pt)
}
