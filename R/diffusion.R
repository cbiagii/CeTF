#' @title Network diffusion analysis
#'
#' @description Expand node selection using network propagation algorithms 
#' generating the expanded network for a core of genes and the 
#' network plot of this subnetwork.
#'
#' @param object CeTF object resulted from \code{\link{runAnalysis}} function.
#' @param cond Which conditions to be used to perform the diffusion analysis. 
#' The options are: network1 (1th condition) and network2 (2th condition).
#' @param genes A single gene or a vector of characters indicating which genes
#' will be used to perform diffusion analysis.
#' @param cyPath System path of Cytoscape software (see \emph{details} for further 
#' informations).
#' @param name Network output name (default: top_diffusion)
#' @param label If label is TRUE, shows the names of nodes (default: TRUE).
#'
#' @return Returns a list with the plot of the network and a table with the 
#' diffusion network.
#'
#' @details 
#' To perform the diffusion analysis is necessary to install the latest Cytoscape 
#' software version (\url{https://cytoscape.org/}).
#' 
#' The \emph{cyPath} argument varies depending on the operating system used, for
#' example:
#' \enumerate{
#'     \item \strong{For Windows users:} C:/Program Files/Cytoscape_v3.7.2/Cytoscape.exe
#'     \item \strong{For Linux users:} /home/user/Cytoscape_v3.7.2/Cytoscape
#'     \item \strong{For macOS users:} 
#' }
#' 
#' @importFrom RCy3 cytoscapePing createNetworkFromIgraph selectNodes commandsPOST createSubnetwork createIgraphFromNetwork
#' @importFrom igraph graph_from_data_frame as_long_data_frame
#' @importFrom network network
#' @importFrom GGally ggnet2  
#' @importFrom ggplot2 coord_equal guides ggtitle
#'
#' @examples
#' \dontrun{ 
#' data(CeTFdemo)
#'
#' result <- diffusion(object = CeTFdemo, 
#'                     cond = 'network1', 
#'                     genes = c('ENSG00000185591', 'ENSG00000179094'), 
#'                     cyPath = 'C:/Program Files/Cytoscape_v3.7.2/Cytoscape.exe', 
#'                     name = 'top_diffusion',
#'                     label = TRUE)
#' }
#'
#' @export
diffusion <- function(object, cond, genes, cyPath, name = "top_diffusion", 
    label = TRUE) {
    if (!is(object, "CeTF")) {
        stop("the input must be a CeTF class object")
    }
    if (missing(cond)) {
        stop("cond argument must be network1 or network2")
    }
    if (missing(genes)) {
        stop("You must input a single gene/TF or a character with genes/TFs")
    }
    if (!file.exists(cyPath)) {
        stop("The path to Cytoscape software must be valid")
    }
    
    if (tryCatch(cytoscapePing(), error = function(e) {
        "wait"
    }) != "You are connected to Cytoscape!") {
        # Openign Cytoscape software
        message("Opening Cytoscape...")
        system2(cyPath, wait = FALSE)
        while (TRUE) {
            if (invisible(tryCatch(cytoscapePing(), error = function(e) {
                "wait"
            }) == "wait")) {
                Sys.sleep(2)
            } else {
                (break)()
            }
        }
        message("Cytoscape opened!!")
    }
    
    ig <- graph_from_data_frame(NetworkData(object, cond), directed = FALSE)
    createNetworkFromIgraph(ig, "network")
    
    selectNodes(nodes = genes, by.col = "name")
    commandsPOST("diffusion diffuse")
    
    createSubnetwork("selected", subnetwork.name = name)
    
    ig2 <- createIgraphFromNetwork(name)
    net <- as_long_data_frame(ig2)[, c("from_name", "to_name")]
    colnames(net) <- c("gene1", "gene2")
    
    # Network plot
    nt1 <- network(net, directed = FALSE)
    genes1 <- unique(c(as.character(net[["gene1"]]), as.character(net[["gene2"]])))
    gns1 <- setdiff(genes1, genes)
    
    map <- data.frame(pattern = c(genes, gns1), sb = c(rep("Selected", 
        length(genes)), rep("Gene", length(gns1))), stringsAsFactors = FALSE)
    
    og <- data.frame(singleValue = network.vertex.names(nt1), stringsAsFactors = FALSE)
    cvt <- vapply(map$pattern, grepl, logical(nrow(og)), og$singleValue)
    og$Category <- map$sb[max.col(cvt, "last")]
    x1 <- factor(og$Category)
    
    nt1 %v% "color" <- as.character(x1)
    y <- c("#4DAF4A", "#E41A1C")
    names(y) = levels(x1)
    
    if (label) {
        label <- genes1
    } else {
        label <- c()
    }
    
    pt <- ggnet2(nt1, color = "color", color.legend = "", palette = y, 
        edge.size = 0.5, edge.color = "gray70", label.size = 1, alpha = 0.75, 
        size = "degree", edge.alpha = 0.5, label = label, legend.position = "bottom") + 
        coord_equal() + guides(size = FALSE) + ggtitle("Diffusion Network")
    
    return(list(plot = pt, network = net))
}
