#' @title Heatmap-like functional classification
#'
#' @description Heatmap-like functional classification to visualize the 
#' enrichment analysis results from \code{\link{getEnrich}} function. The plot
#' contains the heatmap with the associated pathways genes, the significance
#' of the enrichment and a barplot with the enrichment ratio.
#'
#' @param res A dataframe with \code{\link{getEnrich}} results.
#' @param diff A dataframe with all differentialy expressed genes obtained
#' from \code{\link{runAnalysis}} function. For better understanding, simply
#' use the \code{\link{getDE}} accessor with 'all' option.
#' @param showCategory Number of enriched terms to display (default: 10).
#' @param font_size Size of gene row names (default: 6).
#'
#' @return Returns a Heatmap-like functional classification
#'
#' @importFrom utils head
#' @importFrom circlize colorRamp2
#' @importFrom stats dist hclust
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation rowAnnotation Legend 
#' draw max_text_width anno_text anno_barplot anno_simple
#' @importFrom grid gpar unit
#'
#' @examples
#' # loading enrichdemo and CeTFdemo object
#' data(enrichdemo)
#' data(CeTFdemo)
#' 
#' heatPlot(res = enrichdemo$results, 
#'          diff = getDE(CeTFdemo, 'all'), 
#'          showCategory = 10)
#'
#' @export
heatPlot <- function(res, diff, showCategory = 10, font_size = 6) {
    if (missing(res)) {
        stop("res must be a dataframe with enrichment results")
    }
    if (missing(diff)) {
        stop("res must be a dataframe with differential expression values")
    }
    if (showCategory == 1) {
        stop("showCategory must be greater than 1")
    }
    
    res <- head(res, showCategory)
    
    pathways <- res[["ID"]]
    genes <- sort(unique(unlist(strsplit(res[["geneID"]], "\\/"))))
    
    mat <- matrix(0, nrow = length(pathways), ncol = length(genes))
    rownames(mat) <- pathways
    colnames(mat) <- genes
    
    for (i in seq_along(pathways)) {
        gns <- sort(unique(unlist(strsplit(res[i, "geneID"], "\\/"))))
        mat[i, which(colnames(mat) %in% gns)] <- diff[gns, 3]
    }
    
    is_sig = res[["p.adjust"]] < 0.05
    pch = rep("*", length(res[["p.adjust"]]))
    pch[!is_sig] <- NA
    pvalue_col_fun <- colorRamp2(c(0, 1, 2), c("green", "white", "red"))
    
    clust1 <- hclust(dist(mat))
    clust2 <- hclust(dist(t(mat)))
    
    mat[mat == 0] <- NA
    cn = colnames(mat)
    
    breaks <- seq(-4, 4, by = 0.1)
    ht <- Heatmap(mat, show_column_names = FALSE, bottom_annotation = HeatmapAnnotation(text = anno_text(cn, 
        rot = 45, location = unit(1, "npc"), just = "right", gp = gpar(fontsize = font_size)), 
        annotation_height = max_text_width(cn)), cluster_rows = clust1, 
        cluster_columns = clust2, rect_gp = gpar(col = "darkgrey"), heatmap_legend_param = list(nrow = 1, 
            title = "Diff"), show_column_dend = FALSE, show_row_dend = FALSE, 
        row_names_gp = gpar(fontsize = 10))
    
    ann1 <- rowAnnotation(enrichmentRatio = anno_barplot(res[['Count']]/length(genes), 
        width = unit(3, "cm"), axis_param = list(direction = "reverse")))
    ann2 <- rowAnnotation(FDR = anno_simple(-log10(res[["p.adjust"]]), col = pvalue_col_fun, 
        pch = pch))
    
    total <- ann1 + ann2 + ht
    
    lgd_pvalue = Legend(title = "FDR", col_fun = pvalue_col_fun, at = c(0, 
        1, 2), labels = c("0.1", "0.01", "0.001"), title_gp = gpar(fontsize = 8))
    lgd_sig = Legend(pch = "*", type = "points", labels = "< 0.05", title_gp = gpar(fontsize = 8))
    
    return(draw(total, annotation_legend_list = list(lgd_pvalue, lgd_sig), 
        merge_legends = TRUE, heatmap_legend_side = "right"))
}