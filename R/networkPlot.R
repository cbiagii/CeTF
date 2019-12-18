#' @title Network integrated plot
#'
#' @description Generate the plot for network that integrates the network data
#' from a condition with \code{\link{getGroupGO}} results.
#'
#' @param netCond Network of a specific condition. Can be found in result of \code{\link{runAnalysis}} (see \code{\link{getNet1}} and \code{\link{getNet2}})
#' @param netGO Dataframe with the results of \code{\link{getGroupGO}}
#' (second element of list).
#' @param keyTFs TFs identified as importants by \code{\link{runAnalysis}} (see \code{\link{getKeyTF}})
#'
#' @return Returns an network plot with interaction between genes, TFs and/or ontologies.
#'
#' @importFrom network network network.vertex.names '%v%' '%v%<-'
#' @importFrom GGally ggnet2
#' @importFrom utils head tail
#' @importFrom ggplot2 coord_equal guides
#'
#' @examples
#' # load the annotation package
#' library(org.Hs.eg.db)
#'
#' # load the pcitrif class object resulted from runAnalysis function
#' data(pcitrifExample)
#'
#' # getting the genes in network of condition 1
#' genes <- unique(c(as.character(getNet1(pcitrifExample)[,1]),
#'                  as.character(getNet1(pcitrifExample)[,2])))
#'
#' # performing getGroupGO analysis
#' cond1 <- getGroupGO(genes = genes,
#'                     ont = 'BP',
#'                     keyType = 'ENSEMBL',
#'                     annoPkg = org.Hs.eg.db)
#'
#' # selecting only first 12 pathways
#' t1 <- head(cond1$results, 12)
#'
#' # subsetting the network to have only the first 12 pathways
#' t2 <- subset(cond1$netGO, cond1$netGO$gene1 %in% as.character(t1[,1]))
#'
#' # generate the networkPlot for condition1
#' networkPlot(netCond = getNet1(pcitrifExample),
#'             netGO = t2,
#'             keyTFs = getKeyTF(pcitrifExample))
#'
#'
#'
#' @export
networkPlot <- function(netCond, netGO, keyTFs) {
    if(missing(netCond)){stop("No \"netCond\" parameter provided")}
    if(missing(netGO)){stop("No \"netGO\" parameter provided")}
    if(missing(keyTFs)){stop("No \"keyTFs\" parameter provided")}

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
