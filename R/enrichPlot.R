#' @title Plots to visualize the enrichment analysis results
#'
#' @description Generate three types of plots to visualize the enrichment 
#' analysis results from \code{\link{getEnrich}} function. The plots are an 
#' circular barplot, barplot and dotplot.
#'
#' @param res A dataframe with \code{\link{getEnrich}} results.
#' @param showCategory Number of enriched terms to display (default: 10).
#' @param type Type of plot: circle, bar or dot (default: circle).
#'
#' @return Returns a circle, bar or dot plot of enrichment analysis results.
#'
#' @importFrom utils head 
#' @importFrom ggplot2 ggplot aes geom_bar coord_polar theme_minimal theme 
#' element_blank labs coord_flip theme_bw scale_x_discrete scale_fill_gradient 
#' geom_point scale_color_gradient scale_size_continuous scale_fill_continuous
#'
#' @examples
#' # loading enrichdemo
#' data(enrichdemo)
#' 
#' # circle barplot
#' enrichPlot(res = enrichdemo$results, 
#'            showCategory = 10, 
#'            type = 'circle')
#' 
#' # barplot
#' enrichPlot(res = enrichdemo$results, 
#'            showCategory = 10, 
#'            type = 'bar')
#' 
#' # dotplot
#' enrichPlot(res = enrichdemo$results, 
#'            showCategory = 10, 
#'            type = 'dot')
#'
#' @export
enrichPlot <- function(res, showCategory = 10, type = "circle") {
    if (missing(res)) {
        stop("res must be a dataframe with enrichment results")
    }
    if (showCategory == 1) {
        stop("showCategory must be greater than 1")
    }
    
    res <- head(res, showCategory)
    res$enrichmentRatio <- res[['Count']]/length(genes)
    
    if (type == "circle") {
        pt <- ggplot(res, aes(x = .data[["Description"]], y = .data[["enrichmentRatio"]]), 
            colour = .data[["enrichmentRatio"]]) + geom_bar(aes(fill = .data[["enrichmentRatio"]], 
            alpha = .data[["p.adjust"]]), stat = "identity", position = "stack") + 
            coord_polar(clip = "off") + theme_minimal() + theme(axis.title.x = element_blank(), 
            axis.title.y = element_blank(), axis.text.y = element_blank(), 
            axis.text.x = element_text(color = "grey20", size = 8, hjust = 0.5, 
                vjust = 0.5, face = "plain"), legend.position = "bottom") + 
            labs(fill = "enrichmentRatio", alpha = "p.adjust")
    } else if (type == "bar") {
        positions <- res[order(res[["enrichmentRatio"]], decreasing = FALSE), 
            "Description"]
        pt <- ggplot(data = res, aes(x = .data[["Description"]], y = .data[["enrichmentRatio"]])) + 
            geom_bar(stat = "identity", aes(fill = .data[["p.adjust"]])) + coord_flip() + 
            theme_bw() + scale_x_discrete(limits = positions) + scale_fill_gradient(low = "#F4B41A", 
            high = "#143D59", trans = "reverse") + theme(axis.title.x = element_blank(), 
            axis.title.y = element_blank()) + labs(fill = "p.adjust")
    } else if (type == "dot") {
        positions <- res[order(res[["enrichmentRatio"]], decreasing = FALSE), 
            "Description"]
        pt <- ggplot(res, aes(x = .data[["Description"]], y = .data[["enrichmentRatio"]], 
            size = .data[["Count"]], color = .data[["p.adjust"]])) + geom_point(alpha = 0.8) + 
            theme_bw() + scale_color_gradient(low = "#F4B41A", high = "#143D59", 
            trans = "reverse") + scale_x_discrete(limits = positions) + 
            coord_flip() + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + 
            labs(color = "p.adjust", size = "Count") + scale_size_continuous(range = range(res[["Count"]]))
    }
    return(pt)
}
