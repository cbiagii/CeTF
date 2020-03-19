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
#' @param size Size of nodes labels (default: 0.5).
#' @param type Type of plot selected (GO or Integrated). If GO will plot the 
#' associated GO grouped by some variable, and if Integrated will plot a 
#' integrated network with genes, GO and TFs.
#'
#' @return Returns a list with the plot of the network for GO or integrated 
#' output under a condition and the table used to plot the network.
#'
#' @importFrom ggnetwork ggnetwork geom_edges geom_nodes geom_nodetext theme_facet
#' @importFrom ggplot2 fortify ggplot facet_wrap element_rect coord_equal guides scale_color_brewer
#' @importFrom network network network.vertex.names '%v%' '%v%<-'
#' @importFrom GGally ggnet2
#' @importFrom utils head tail
#' @importFrom network set.edge.attribute
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
#'               keyTFs = NetworkData(CeTFdemo, 'keytfs'), 
#'               type = 'Integrated')
#' pt$plot
#' head(pt$tab)
#' }
#'
#' @export
netGOTFPlot <- function(netCond, resultsGO, netGO, anno, groupBy = "pathways", 
    TFs = NULL, genes = NULL, keyTFs = NULL, size = 0.5, type = NULL) {
    
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
                gns <- as.character(subset(netGO, netGO[["gene1"]] == as.character(x[["ID"]]))[, 
                  "gene2"])
                tmp1 <- netCond[which(netCond[["gene1"]] %in% gns & netCond[["gene2"]] %in% 
                  gns), ]
                if (nrow(tmp1) != 0) {
                  tmp1[["pathway"]] <- as.character(x[["ID"]])
                }
                return(tmp1)
            })
            tmp <- do.call(rbind, tmp)
            
            # adapted fortify function
            allnodes <- expand.grid(unique(anno[["genes"]]), unique(tmp[["pathway"]]), 
                stringsAsFactors = FALSE)
            node.data.expanded <- merge(allnodes, anno, by.x = "Var1", 
                by.y = "genes")
            tab <- merge(tmp, node.data.expanded, by.x = c("gene1", "pathway"), 
                by.y = c("Var1", "Var2"), all = TRUE)
            
            tab1 <- na.omit(tab)
            tab1 <- unique(tab1)
            
            em.path <- lapply(unique(tab1[["pathway"]]), function(x) subset(tab1, 
                pathway == x)[, c(1, 3)])
            em.pathway <- lapply(em.path, network, directed = FALSE)
            em.tmp <- lapply(unique(tab[["pathway"]]), function(x) subset(tab, 
                pathway == x))
            
            for (i in seq_along(em.pathway)) {
                set.edge.attribute(em.pathway[[i]], "pathway", em.tmp[[i]][["pathway"]])
                
                df <- unique(data.frame(genes = as.character(em.tmp[[i]][["gene1"]]), 
                  class = as.character(em.tmp[[i]][["class"]]), stringsAsFactors = FALSE))
                var1 <- df[["class"]]
                names(var1) <- df[["genes"]]
                em.pathway[[i]] %v% "class" <- var1[network.vertex.names(em.pathway[[i]])]
            }
            
            g <- list()
            for (i in seq_along(em.pathway)) {
                g[[i]] <- ggplot(ggnetwork(em.pathway[[i]], arrow.gap = 0.02, 
                  by = "pathway", layout = "fruchtermanreingold"), aes(x, 
                  y, xend = xend, yend = yend)) + geom_edges(aes(color = class), 
                  alpha = 0.2, curvature = 0.05, color = "gray55") + geom_nodes(aes(color = class), 
                  size = 3, show.legend = TRUE) + geom_nodetext(aes(label = vertex.names), 
                  size = size) + scale_color_brewer("Type", palette = "Set1") + 
                  facet_wrap(~pathway) + theme_facet(legend.position = "bottom")
            }
            
            pt <- ggarrange(plotlist = g, common.legend = TRUE, legend = "bottom")
            
            out <- list(plot = pt, tab = suppressWarnings(split(tab1, f = tab[["pathway"]])))
        } else if (groupBy == "TFs") {
            tmp1 <- apply(resultsGO, 1, function(x) {
                gns <- as.character(subset(netGO, netGO[["gene1"]] == as.character(x[["ID"]]))[, 
                  "gene2"])
                tmp1 <- data.frame(genes = gns)
                tmp1[["pathway"]] <- as.character(x[["ID"]])
                return(tmp1)
            })
            tmp1 <- do.call(rbind, tmp1)
            
            if (length(which(TFs %in% as.character(tmp1[["genes"]]))) == 
                0) {
                stop("there isn't none TF passed as input in the groupGO analysis: select others TFs")
            }
            
            tmp2 <- NULL
            for (i in seq_len(length(TFs))) {
                t1 <- subset(tmp1, tmp1[["genes"]] == TFs[i])
                if (nrow(t1) != 0) {
                  t2 <- subset(netGO, netGO[["gene1"]] %in% t1[["pathway"]])
                  t2[["TF"]] <- TFs[i]
                  tmp2 <- rbind(tmp2, t2)
                } else {
                  next
                }
            }
            
            anno <- rbind(anno, data.frame(genes = unique(as.character(tmp2[["gene1"]])), 
                class = "pathway"))
            
            # adapted fortify function
            allnodes <- expand.grid(unique(anno[["genes"]]), unique(tmp2[["TF"]]), 
                stringsAsFactors = FALSE)
            node.data.expanded <- merge(allnodes, anno, by.x = "Var1", 
                by.y = "genes")
            tab <- merge(tmp2, node.data.expanded, by.x = c("gene1", "TF"), 
                by.y = c("Var1", "Var2"), all = TRUE)
            
            tab1 <- na.omit(tab)
            tab1 <- unique(tab1)
            
            em.1 <- lapply(unique(tab1[["TF"]]), function(x) subset(tab1, 
                TF == x)[, c(1, 3)])
            em.TF <- lapply(em.1, network, directed = FALSE)
            em.tmp <- lapply(unique(tab[["TF"]]), function(x) subset(tab, 
                TF == x))
            
            for (i in seq_along(em.TF)) {
                set.edge.attribute(em.TF[[i]], "TF", em.tmp[[i]][["TF"]])
                
                df <- unique(data.frame(genes = as.character(em.tmp[[i]][["gene1"]]), 
                  class = as.character(em.tmp[[i]][["class"]]), stringsAsFactors = FALSE))
                var1 <- df[["class"]]
                names(var1) <- df[["genes"]]
                em.TF[[i]] %v% "class" <- var1[network.vertex.names(em.TF[[i]])]
            }
            
            g <- list()
            for (i in seq_along(em.TF)) {
                g[[i]] <- ggplot(ggnetwork(em.TF[[i]], arrow.gap = 0.02, 
                  by = "TF", layout = "fruchtermanreingold"), aes(x, y, 
                  xend = xend, yend = yend)) + geom_edges(aes(color = class), 
                  alpha = 0.2, curvature = 0.05, color = "gray55") + geom_nodes(aes(color = class), 
                  size = 3, show.legend = TRUE) + geom_nodetext(aes(label = vertex.names), 
                  size = size) + scale_color_brewer("Type", palette = "Set1") + 
                  facet_wrap(~TF) + theme_facet(legend.position = "bottom")
            }
            
            pt <- ggarrange(plotlist = g, common.legend = TRUE, legend = "bottom")
            
            out <- list(plot = pt, tab = suppressWarnings(split(tab, f = tab[["TF"]])))
        } else if (groupBy == "genes") {
            tmp1 <- apply(resultsGO, 1, function(x) {
                gns <- as.character(subset(netGO, netGO[["gene1"]] == as.character(x[["ID"]]))[, 
                  "gene2"])
                tmp1 <- data.frame(genes = gns)
                tmp1[["pathway"]] <- as.character(x[["ID"]])
                return(tmp1)
            })
            tmp1 <- do.call(rbind, tmp1)
            
            tmp2 <- NULL
            for (i in seq_len(length(genes))) {
                t1 <- subset(tmp1, tmp1[["genes"]] == genes[i])
                t2 <- subset(netGO, netGO[["gene1"]] %in% t1[["pathway"]])
                t2[["genes"]] <- genes[i]
                tmp2 <- rbind(tmp2, t2)
            }
            
            anno <- rbind(anno, data.frame(genes = unique(as.character(tmp2[["gene1"]])), 
                class = "pathway"))
            
            # adapted fortify function
            allnodes <- expand.grid(unique(anno[["genes"]]), unique(tmp2[["genes"]]), 
                stringsAsFactors = FALSE)
            node.data.expanded <- merge(allnodes, anno, by.x = "Var1", 
                by.y = "genes")
            tab <- merge(tmp2, node.data.expanded, by.x = c("gene1", "genes"), 
                by.y = c("Var1", "Var2"), all = TRUE)
            
            tab1 <- na.omit(tab)
            tab1 <- unique(tab1)
            
            em.1 <- lapply(unique(tab1[["genes"]]), function(x) subset(tab1, 
                genes == x)[, c(1, 3)])
            em.genes <- lapply(em.1, network, directed = FALSE)
            em.tmp <- lapply(unique(tab[["genes"]]), function(x) subset(tab, 
                genes == x))
            
            for (i in seq_along(em.genes)) {
                set.edge.attribute(em.genes[[i]], "genes", em.tmp[[i]][["genes"]])
                
                df <- unique(data.frame(genes = as.character(em.tmp[[i]][["gene1"]]), 
                  class = as.character(em.tmp[[i]][["class"]]), stringsAsFactors = FALSE))
                var1 <- df[["class"]]
                names(var1) <- df[["genes"]]
                em.genes[[i]] %v% "class" <- var1[network.vertex.names(em.genes[[i]])]
            }
            
            g <- list()
            for (i in seq_along(em.genes)) {
                g[[i]] <- ggplot(ggnetwork(em.genes[[i]], arrow.gap = 0.02, 
                  by = "genes", layout = "fruchtermanreingold"), aes(x, 
                  y, xend = xend, yend = yend)) + geom_edges(aes(color = class), 
                  alpha = 0.2, curvature = 0.05, color = "gray55") + geom_nodes(aes(color = class), 
                  size = 3, show.legend = TRUE) + geom_nodetext(aes(label = vertex.names), 
                  size = size) + scale_color_brewer("Type", palette = "Set1") + 
                  facet_wrap(~genes) + theme_facet(legend.position = "bottom")
            }
            
            pt <- ggarrange(plotlist = g, common.legend = TRUE, legend = "bottom")
            
            out <- list(plot = pt, tab = suppressWarnings(split(tab, f = tab[["genes"]])))
        }
    } else if (type == "Integrated") {
        network <- rbind(netCond, netGO)
        net <- network(network, directed = FALSE)
        
        values <- unique(c(as.character(network[["gene1"]]), as.character(network[["gene2"]])))
        pathways <- unique(as.character(netGO[["gene1"]]))
        
        if (is.character(keyTFs)) {
            TFs <- as.character(keyTFs)
        } else {
            TFs <- as.character(keyTFs[["TF"]])
        }
        gns <- setdiff(values, c(TFs, pathways))
        
        map <- data.frame(pattern = c(TFs, gns, pathways), sb = c(rep("TFs", 
            length(TFs)), rep("Gene", length(gns)), rep("GO", length(pathways))), 
            stringsAsFactors = FALSE)
        og <- data.frame(singleValue = network.vertex.names(net), stringsAsFactors = FALSE)
        cvt <- vapply(map[["pattern"]], grepl, logical(nrow(og)), og[["singleValue"]])
        og[["Category"]] <- map[["sb"]][max.col(cvt, "last")]
        x <- factor(og[["Category"]])
        
        net %v% "color" <- net %v% "class" <- as.character(x)
        
        pt <- ggplot(ggnetwork(net, by = "color", layout = "spring"), aes(x, 
            y, xend = xend, yend = yend)) + geom_edges(aes(color = class), 
            alpha = 0.2, curvature = 0.05, color = "gray55") + geom_nodes(aes(color = class), 
            size = 3, show.legend = TRUE) + geom_nodetext(aes(label = vertex.names), 
            size = size) + scale_color_brewer("Type", palette = "Set1") + 
            theme_facet(legend.position = "bottom")
        
        out <- list(plot = pt, tab = network)
    }
    return(out)
}
