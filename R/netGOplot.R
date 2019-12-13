#' @title getGroupGO network plot
#'
#' @description Generate the plot of groupGO network result.
#'
#' @param netCond Network of a specific condition. Can be found in result of \link[pcitRif]{runAnalysis} (step4 -> network_cond1 or network_cond2)
#' @param resultsGO Dataframe with the results of \link[pcitRif]{getGroupGO} (first element of list). This result can be decreased applying filters for the pathways selection.
#' @param netGO Dataframe with the results of \link[pcitRif]{getGroupGO} (second element of list).
#' @param anno Annotation of gene or TFs. Can be found in result of runAnalysis (step4 -> anno)
#' @param label If label is TRUE, shows the names of nodes.
#'
#' @return The network for groupGO for a determined condition.
#'
#' @importFrom geomnet geom_net theme_net as.edgedf
#' @importFrom ggplot2 fortify ggplot facet_wrap
#'
#' @examples
#' \dontrun{
#' #See vignette for more details
#' netGOplot(netCond,
#' resultsGO,
#' netGO,
#' anno,
#' label = F)
#' }
#'
#'
#'
#' @export
netGOplot <- function(netCond, resultsGO, netGO, anno, label = F) {
  tmp <- NULL
  for (i in 1:nrow(resultsGO)) {
    genes <- as.character(subset(netGO, netGO$gene1 == as.character(resultsGO$ID[i]))[,2])
    tmp1 <- netCond[which(netCond$gene1 %in% genes & netCond$gene2 %in% genes), ]

    if (nrow(tmp1) != 0) {
      tmp1$pathway <- as.character(resultsGO$ID[i])
      tmp <- rbind(tmp, tmp1)
    } else { next }
  }

  data <- fortify(as.edgedf(tmp), anno, group = "pathway")
  pt <- ggplot(data, aes(from_id = from, to_id = to_id)) +
    geom_net(aes(colour = class, group = class, linewidth = 3 * (...samegroup.. / 8 + .125)),
             layout.alg = "fruchtermanreingold",
             ealpha = 0.5, size = 3, curvature = 0.05,
             directed = F, arrowsize = 0.5, show.legend = T,
             fiteach = T, labelon = label, fontsize=0.5, alpha = 0.25,
             labelcolour = "black", singletons = FALSE) +
    facet_wrap(~pathway) +
    theme_net() + theme(panel.background = element_rect(colour = 'black'))

  return(pt)
}
