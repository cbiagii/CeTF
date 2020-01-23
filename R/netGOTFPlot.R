#' @title Plot a network for Ontologies, genes and TFs
#'
#' @description Generate the plot of groupGO network result of
#' \code{\link{getGroupGO}} function, and the integrated network of genes, GOs 
#' and TFs.
#'
#' @param netCond Network of a specific condition. Can be found in
#' result of \code{\link{runAnalysis}} (see \code{\link{NetworkData}} and
#' \code{\link{NetworkData}}).
#' @param resultsGO Dataframe with the results of \code{\link{getGroupGO}}
#' (first element of list). This result can be filtered by applying filters
#' for pathways selection.
#' @param netGO Dataframe with the results of \code{\link{getGroupGO}}
#' (second element of list).
#' @param anno Annotation of gene or TFs. Can be found in result of
#' \code{\link{runAnalysis}} function (see \code{\link{NetworkData}}).
#' @param groupBy Which variables do you want to group in GO type? The options 
#' are: 'pathways', 'TFs' and 'genes' (default: 'pathways').
#' @param TFs A character with selected TFs.
#' @param genes A character with selected genes.
#' @param keyTFs TFs identified as importants by \code{\link{runAnalysis}} 
#' (see \code{\link{NetworkData}}). This argument is used only if the type argument
#' equals Integrated.
#' @param label If label is TRUE, shows the names of nodes (default: FALSE).
#' @param type Type of plot selected (GO or Integrated). If GO will plot the 
#' associated GO grouped by some variable, and if Integrated will plot a 
#' integrated network with genes, GO and TFs.
#'
#' @return Returns a list with the plot of the network for GO or integrated 
#' output under a condition and the table used to plot the network.
#'
#' @importFrom geomnet geom_net theme_net as.edgedf
#' @importFrom ggplot2 fortify ggplot facet_wrap element_rect coord_equal guides
#' @importFrom network network network.vertex.names '%v%' '%v%<-'
#' @importFrom GGally ggnet2
#' @importFrom utils head tail
#'
#' @examples
#' \dontrun{ 
#' # load the annotation package
#' library(org.Hs.eg.db)
#'
#' # load the CeTF class object resulted from runAnalysis function
#' data(CeTFdemo)
#'
#' # getting the genes in network of condition 1
#' genes <- unique(c(as.character(NetworkData(CeTFdemo, 'network1')[, 'gene1']),
#'                  as.character(NetworkData(CeTFdemo, 'network1')[, 'gene2'])))
#'
#' # performing getGroupGO analysis
#' cond1 <- getGroupGO(genes = genes,
#'                     ont = 'BP',
#'                     keyType = 'ENSEMBL',
#'                     annoPkg = org.Hs.eg.db, 
#'                     level = 3)
#'
#' # selecting only first 12 pathways
#' t1 <- head(cond1$results, 12)
#'
#' # subsetting the network to have only the first 12 pathways
#' t2 <- subset(cond1$netGO, cond1$netGO$gene1 %in% as.character(t1[, 'ID']))
#'
#' # generate the GO plot grouping by pathways
#' pt <- netGOTFPlot(netCond = NetworkData(CeTFdemo, 'network1'),
#'               resultsGO = t1,
#'               netGO = t2,
#'               anno = NetworkData(CeTFdemo, 'annotation'),
#'               groupBy = 'pathways',
#'               label = TRUE, 
#'               keyTFs = NetworkData(CeTFdemo, 'keytfs'), 
#'               type = 'GO')
#' pt$plot
#' head(pt$tab$`GO:0006807`)
#' 
#' # generate the Integrated plot
#' pt <- netGOTFPlot(netCond = NetworkData(CeTFdemo, 'network1'),
#'               resultsGO = t1,
#'               netGO = t2,
#'               anno = NetworkData(CeTFdemo, 'annotation'),
#'               groupBy = 'pathways',
#'               label = TRUE, 
#'               keyTFs = NetworkData(CeTFdemo, 'keytfs'), 
#'               type = 'Integrated')
#' pt$plot
#' head(pt$tab)
#' }
#'
#' @export
netGOTFPlot <- function(netCond, resultsGO, netGO, anno, groupBy = "pathways", 
    TFs = NULL, genes = NULL, keyTFs = NULL, label = FALSE, type = NULL) {
    if (groupBy == "TFs" & is.null(TFs)) {
        stop("for groupBy = TFs you must input some TFs")
    }
    if (groupBy == "genes" & is.null(genes)) {
        stop("for groupBy = genes you must input some genes")
    }
    if (type == "Integrated" & is.null(keyTFs)) {
        stop("this argument must be a character of ketTFs outputed from runAnalysis function")
    }
    if (is.null(type)) {
        stop("this argument must be one of options: GO or all")
    }
    
    if (type == "GO") {
        if (groupBy == "pathways") {
            tmp <- apply(resultsGO, 1, function(x) {
                gns <- as.character(subset(netGO, netGO$gene1 == as.character(x[["ID"]]))[, 
                  "gene2"])
                tmp1 <- netCond[which(netCond$gene1 %in% gns & netCond$gene2 %in% 
                  gns), ]
                if (nrow(tmp1) != 0) {
                    tmp1$pathway <- as.character(x[["ID"]])
                }
                return(tmp1)
            })
            tmp <- do.call(rbind, tmp)
            
            suppressMessages(tab <- fortify(as.edgedf(tmp), anno, group = "pathway"))
            pt <- ggplot(tab, aes(from_id = tab[["from"]], to_id = tab[["to_id"]])) + 
                geom_net(aes(colour = class, group = class, linewidth = 0.5), 
                         layout.alg = "fruchtermanreingold", ealpha = 0.5, 
                         size = 3, curvature = 0.05, directed = FALSE, 
                         arrowsize = 0.5, na.rm = TRUE, show.legend = TRUE, 
                         fiteach = TRUE, labelon = label, fontsize = 0.5, 
                         alpha = 0.25, labelcolour = "black", singletons = FALSE) + 
                theme_net()
            
            if (length(unique(tab[["pathway"]])) == 1) {
                pt <- pt + ggtitle(paste("Network for", unique(tab[["pathway"]])))
            } else {
                pt <- pt + facet_wrap(~pathway) + 
                    theme(panel.background = element_rect(colour = "black")) + 
                    ggtitle("Network for Pathways")
            }
            out <- list(plot = pt, tab = split(tab, f = tab$pathway))
        } else if (groupBy == "TFs") {
            tmp1 <- apply(resultsGO, 1, function(x) {
                gns <- as.character(subset(netGO, netGO$gene1 == as.character(x[["ID"]]))[, 
                  "gene2"])
                tmp1 <- data.frame(genes = gns)
                tmp1$pathway <- as.character(x[["ID"]])
                return(tmp1)
            })
            tmp1 <- do.call(rbind, tmp1)
            
            if (length(which(TFs %in% as.character(tmp1$genes))) == 0) {
                stop("there isn't none TF passed as input in the groupGO analysis: select others TFs")
            }
            
            tmp2 <- NULL
            for (i in seq_len(length(TFs))) {
                t1 <- subset(tmp1, tmp1$genes == TFs[i])
                if (nrow(t1) != 0) {
                  t2 <- subset(netGO, netGO$gene1 %in% t1$pathway)
                  t2$TF <- TFs[i]
                  tmp2 <- rbind(tmp2, t2)
                } else {
                  next
                }
            }
            
            anno <- rbind(anno, data.frame(genes = unique(as.character(tmp2$gene1)), 
                class = "pathway"))
            suppressMessages(tab <- fortify(as.edgedf(tmp2), anno, group = "TF"))
            pt <- ggplot(tab, aes(from_id = tab[["from"]], to_id = tab[["to_id"]])) + 
                geom_net(aes(colour = class, group = class, linewidth = 0.5), 
                         layout.alg = "fruchtermanreingold", ealpha = 0.5, size = 3, 
                         curvature = 0.05, directed = FALSE, arrowsize = 0.5, 
                         na.rm = TRUE, show.legend = TRUE, fiteach = TRUE, labelon = label, 
                         fontsize = 0.5, alpha = 0.25, labelcolour = "black", 
                         singletons = FALSE) + theme_net()
            
            if (length(unique(tab[["TF"]])) == 1) {
                pt <- pt + ggtitle(paste("Network for", unique(tab[["TF"]])))
            } else {
                pt <- pt + facet_wrap(~TF) + 
                    theme(panel.background = element_rect(colour = "black")) + 
                    ggtitle("Network for Pathways")
            }
            out <- list(plot = pt, tab = split(tab, f = tab$TF))
        } else if (groupBy == "genes") {
            tmp1 <- apply(resultsGO, 1, function(x) {
                gns <- as.character(subset(netGO, netGO$gene1 == as.character(x[["ID"]]))[, 
                  "gene2"])
                tmp1 <- data.frame(genes = gns)
                tmp1$pathway <- as.character(x[["ID"]])
                return(tmp1)
            })
            tmp1 <- do.call(rbind, tmp1)
            
            tmp2 <- NULL
            for (i in seq_len(length(genes))) {
                t1 <- subset(tmp1, tmp1$genes == genes[i])
                t2 <- subset(netGO, netGO$gene1 %in% t1$pathway)
                t2$genes <- genes[i]
                tmp2 <- rbind(tmp2, t2)
            }
            
            anno <- rbind(anno, data.frame(genes = unique(as.character(tmp2$gene1)), 
                class = "pathway"))
            suppressMessages(tab <- fortify(as.edgedf(tmp2), anno, group = "genes"))
            pt <- ggplot(tab, aes(from_id = tab[["from"]], to_id = tab[["to_id"]])) + 
                geom_net(aes(colour = class, group = class, linewidth = 0.5), 
                  layout.alg = "fruchtermanreingold", ealpha = 0.5, size = 3, 
                  curvature = 0.05, directed = FALSE, arrowsize = 0.5, 
                  na.rm = TRUE, show.legend = TRUE, fiteach = TRUE, labelon = label, 
                  fontsize = 0.5, alpha = 0.25, labelcolour = "black", 
                  singletons = FALSE) + theme_net()
            if (length(unique(tab[["genes"]])) == 1) {
                pt <- pt + ggtitle(paste("Network for", unique(tab[["genes"]])))
            } else {
                pt <- pt + facet_wrap(~genes) + 
                    theme(panel.background = element_rect(colour = "black")) + 
                    ggtitle("Network for Genes")
            }
            out <- list(plot = pt, tab = split(tab, f = tab$genes))
        }
    } else if (type == "Integrated") {
        network <- rbind(netCond, netGO)
        net <- network(network, directed = FALSE)
        
        values <- unique(c(as.character(network$gene1), as.character(network$gene2)))
        pathways <- unique(as.character(netGO$gene1))
        
        if (is.character(keyTFs)) {
            TFs <- as.character(keyTFs)
        } else {
            TFs <- as.character(keyTFs$TF)
        }
        gns <- setdiff(values, c(TFs, pathways))
        
        map <- data.frame(pattern = c(TFs, gns, pathways), sb = c(rep("TFs", 
            length(TFs)), rep("Gene", length(gns)), rep("GO", length(pathways))), 
            stringsAsFactors = FALSE)
        og <- data.frame(singleValue = network.vertex.names(net), stringsAsFactors = FALSE)
        cvt <- vapply(map$pattern, grepl, logical(nrow(og)), og$singleValue)
        og$Category <- map$sb[max.col(cvt, "last")]
        x <- factor(og$Category)
        
        net %v% "color" = as.character(x)
        y <- c("#4DAF4A", "#E41A1C", "#377EB8")
        names(y) = levels(x)
        
        pt <- ggnet2(net, color = "color", color.legend = "", palette = y, 
            edge.size = 0.5, edge.color = "gray70", label.size = 1, alpha = 0.75, 
            size = "degree", edge.alpha = 0.5, label = label, 
            legend.position = "bottom", mode = "spring") + coord_equal() + 
            guides(size = FALSE)
        
        out <- list(plot = pt, tab = network)
    }
    return(out)
}
