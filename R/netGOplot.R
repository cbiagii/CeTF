#' @title Plot a network for getGroupGO function
#'
#' @description Generate the plot of groupGO network result of
#' \code{\link{getGroupGO}} function.
#'
#' @param netCond Network of a specific condition. Can be found in
#' result of \code{\link{runAnalysis}} (see \code{\link{getNet1}} and
#' \code{\link{getNet2}}).
#' @param resultsGO Dataframe with the results of \code{\link{getGroupGO}}
#' (first element of list). This result can be filtered by applying filters
#' for pathways selection.
#' @param netGO Dataframe with the results of \code{\link{getGroupGO}}
#' (second element of list).
#' @param anno Annotation of gene or TFs. Can be found in result of
#' \code{\link{runAnalysis}} function (see \code{\link{getAnno}}).
#' @param groupBy Which variables do you want to group? The options are: 'pathways', 'TFs' and 'genes'.
#' @param TFs A character with selected TFs.
#' @param genes A character with selected genes.
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
#' # generate the netGOplot grouping by pathways
#' netGOplot(netCond = getNet1(pcitrifExample),
#'           resultsGO = t1,
#'           netGO = t2,
#'           anno = getAnno(pcitrifExample),
#'           groupBy = 'pathways',
#'           label = TRUE)
#'
#' # generate the netGOplot grouping by TFs
#' TFs <- getKeyTF(pcitrifExample)$TF[1:8]
#' netGOplot(netCond = getNet1(pcitrifExample),
#'           resultsGO = t1,
#'           netGO = t2,
#'           anno = getAnno(pcitrifExample),
#'           groupBy = 'TFs',
#'           TFs = TFs,
#'           label = TRUE)
#'
#' # generate the netGOplot grouping by specifics genes
#' netGOplot(netCond = getNet1(pcitrifExample),
#'           resultsGO = t1,
#'           netGO = t2,
#'           anno = getAnno(pcitrifExample),
#'           groupBy = 'genes',
#'           genes = c('ENSG00000011465', 'ENSG00000026025', 'ENSG00000075624',
#'           'ENSG00000077942', 'ENSG00000087086', 'ENSG00000087245',
#'           'ENSG00000091136', 'ENSG00000100234'),
#'           label = TRUE)
#'
#'
#'
#' @export
netGOplot <- function(netCond, resultsGO, netGO, anno,
    groupBy = "pathways", TFs = NULL, genes = NULL,
    label = FALSE) {
    if (groupBy == "TFs" & is.null(TFs)) {
        stop("for TFs groupBy parameter you must input some TFs")
    }

    if (groupBy == "genes" & is.null(genes)) {
        stop("for genes groupBy parameter you must input some genes")
    }

    if (groupBy == "pathways") {
        tmp <- apply(resultsGO, 1, function(x) {
            gns <- as.character(subset(netGO, netGO$gene1 ==
                as.character(x[["ID"]]))[, 2])
            tmp1 <- netCond[which(netCond$gene1 %in%
                gns & netCond$gene2 %in% gns), ]
            tmp1$pathway <- as.character(x[["ID"]])
            return(tmp1)
        })
        tmp <- do.call(rbind, tmp)

        suppressMessages(tab <- fortify(as.edgedf(tmp),
            anno, group = "pathway"))
        pt <- ggplot(tab, aes(from_id = from, to_id = to_id)) +
            geom_net(aes(colour = class, group = class,
                linewidth = 0.5), layout.alg = "fruchtermanreingold",
                ealpha = 0.5, size = 3, curvature = 0.05,
                directed = FALSE, arrowsize = 0.5, na.rm = TRUE,
                show.legend = TRUE, fiteach = TRUE,
                labelon = label, fontsize = 0.5, alpha = 0.25,
                labelcolour = "black", singletons = FALSE,) +
            ggtitle("Network for Pathways") + facet_wrap(~pathway) +
            theme_net() + theme(panel.background = element_rect(colour = "black"))
    } else if (groupBy == "TFs") {
        tmp1 <- apply(resultsGO, 1, function(x) {
            gns <- as.character(subset(netGO, netGO$gene1 ==
                as.character(x[["ID"]]))[, 2])
            tmp1 <- data.frame(genes = gns)
            tmp1$pathway <- as.character(x[["ID"]])
            return(tmp1)
        })
        tmp1 <- do.call(rbind, tmp1)

        tmp2 <- NULL
        for (i in seq_len(length(TFs))) {
            t1 <- subset(tmp1, tmp1$genes == TFs[i])
            t2 <- subset(netGO, netGO$gene1 %in%
                t1$pathway)
            t2$TF <- TFs[i]
            tmp2 <- rbind(tmp2, t2)
        }

        anno <- rbind(anno, data.frame(genes = unique(as.character(tmp2$gene1)),
            class = "pathway"))
        suppressMessages(tab <- fortify(as.edgedf(tmp2),
            anno, group = "TF"))
        pt <- ggplot(tab, aes(from_id = from, to_id = to_id)) +
            geom_net(aes(colour = class, group = class,
                linewidth = 0.5), layout.alg = "fruchtermanreingold",
                ealpha = 0.5, size = 3, curvature = 0.05,
                directed = FALSE, arrowsize = 0.5, na.rm = TRUE,
                show.legend = TRUE, fiteach = TRUE,
                labelon = label, fontsize = 0.5, alpha = 0.25,
                labelcolour = "black", singletons = FALSE) +
            ggtitle("Network for TFs") + facet_wrap(~TF) +
            theme_net() + theme(panel.background = element_rect(colour = "black"))
    } else if (groupBy == "genes") {
        tmp1 <- apply(resultsGO, 1, function(x) {
            gns <- as.character(subset(netGO, netGO$gene1 ==
                as.character(x[["ID"]]))[, 2])
            tmp1 <- data.frame(genes = gns)
            tmp1$pathway <- as.character(x[["ID"]])
            return(tmp1)
        })
        tmp1 <- do.call(rbind, tmp1)

        tmp2 <- NULL
        for (i in seq_len(length(genes))) {
            t1 <- subset(tmp1, tmp1$genes == genes[i])
            t2 <- subset(netGO, netGO$gene1 %in%
                t1$pathway)
            t2$genes <- genes[i]
            tmp2 <- rbind(tmp2, t2)
        }

        anno <- rbind(anno, data.frame(genes = unique(as.character(tmp2$gene1)),
            class = "pathway"))
        suppressMessages(tab <- fortify(as.edgedf(tmp2),
            anno, group = "genes"))
        pt <- ggplot(tab, aes(from_id = from, to_id = to_id)) +
            geom_net(aes(colour = class, group = class,
                linewidth = 0.5), layout.alg = "fruchtermanreingold",
                ealpha = 0.5, size = 3, curvature = 0.05,
                directed = FALSE, arrowsize = 0.5, na.rm = TRUE,
                show.legend = TRUE, fiteach = TRUE, labelon = label,
                fontsize = 0.5, alpha = 0.25, labelcolour = "black",
                singletons = FALSE) + ggtitle("Network for Genes") +
            facet_wrap(~genes) + theme_net() + theme(panel.background = element_rect(colour = "black"))
    }
    return(pt)
}

