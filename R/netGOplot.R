#' @title Plot a network for getGroupGO function
#'
#' @description Generate the plot of groupGO network result of
#' \code{\link{getGroupGO}} function.
#'
#' @param netCond Network of a specific condition. Can be found in
#' result of \code{\link{runAnalysis}} (see \code{\link{getNet1}} and
#' \code{\link{getNet2}})
#' @param resultsGO Dataframe with the results of \code{\link{getGroupGO}}
#' (first element of list). This result can be filtered by applying filters
#' for pathways selection.
#' @param netGO Dataframe with the results of \code{\link{getGroupGO}}
#' (second element of list).
#' @param anno Annotation of gene or TFs. Can be found in result of
#' \code{\link{runAnalysis}} function (see \code{\link{getAnno}})
#' @param label If label is TRUE, shows the names of nodes.
#'
#' @return The network for \code{\link{getGroupGO}} output under a condition.
#'
#' @importFrom geomnet geom_net theme_net as.edgedf
#' @importFrom ggplot2 fortify ggplot facet_wrap element_rect
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
#'                     ont = "BP",
#'                     keyType = "ENSEMBL",
#'                     annoPkg = org.Hs.eg.db)
#'
#' # selecting only first 12 pathways
#' t1 <- head(cond1$results, 12)
#'
#' # subsetting the network to have only the first 12 pathways
#' t2 <- subset(cond1$netGO, cond1$netGO$gene1 %in% as.character(t1[,1]))
#'
#' # generate the netGOplot
#' netGOplot(netCond = getNet1(pcitrifExample),
#'           resultsGO = t1,
#'           netGO = t2,
#'           anno = getAnno(pcitrifExample),
#'           label = TRUE)
#'
#'
#'
#'
#' @export
netGOplot <- function(netCond, resultsGO, netGO, anno,
    label = FALSE) {
    tmp <- apply(resultsGO, 1, function(x) {
        genes <- as.character(subset(netGO, netGO$gene1 ==
            as.character(x[["ID"]]))[, 2])
        tmp1 <- netCond[which(netCond$gene1 %in% genes &
            netCond$gene2 %in% genes), ]
        tmp1$pathway <- as.character(x[["ID"]])
        return(tmp1)
    })
    tmp <- do.call(rbind, tmp)

    tab <- fortify(as.edgedf(tmp), anno, group = "pathway")
    pt <- ggplot(tab, aes(from_id = from, to_id = to_id)) +
        geom_net(aes(colour = class, group = class,
            linewidth = 0.5), layout.alg = "fruchtermanreingold",
            ealpha = 0.5, size = 3, curvature = 0.05,
            directed = FALSE, arrowsize = 0.5, show.legend = TRUE,
            fiteach = TRUE, labelon = label, fontsize = 0.5,
            alpha = 0.25, labelcolour = "black", singletons = FALSE) +
        facet_wrap(~pathway) + theme_net() + theme(panel.background = element_rect(colour = "black"))

    return(pt)
}
