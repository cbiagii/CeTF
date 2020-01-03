#' @title Histogram of connectivity distribution
#'
#' @description Generate the histogram for adjacency matrix to
#' show the clustering coefficient distribution.
#'
#' @param mat Adjacency matrix resulting from PCIT analysis in which
#' has some zero values.
#'
#' @return Returns the histogram of connectivity distribution.
#'
#' @importFrom ggplot2 ggplot geom_histogram xlab ylab element_line element_blank geom_step aes scale_x_continuous scale_y_continuous element_text
#' @importFrom graphics hist
#' @importFrom ggpubr ggarrange
#' @importFrom scales percent
#'
#' @examples
#' # loading a simulated normalized data
#' data(simNorm)
#'
#' # getting the PCIT results
#' results <- PCIT(simNorm)
#'
#' # plotting the histogram for PCIT significance results
#' histPlot(results$adj_sig)
#'
#'
#'
#' @export
histPlot <- function(mat) {
    if (!is.data.frame(mat) & !is.matrix(mat)) {
        stop("input must be a dataframe or a matrix")
    }
    
    cc <- clustCoef(mat)
    
    df1 <- data.frame(clustcoef = cc)
    pt1 <- ggplot(df1, aes(df1[["clustcoef"]])) + geom_histogram(breaks = seq(0, 
        0.6, by = 0.05), col = "black", fill = "#1F3552") + ggtitle("Connectivity Distribution") + 
        xlab("Proportion of Connections") + ylab("Number of Genes") + scale_y_continuous(expand = c(0, 
        0)) + scale_x_continuous(labels = percent) + theme_bw() + theme(axis.line = element_line(size = 1, 
        colour = "black"), panel.grid.major = element_line(colour = "#d3d3d3"), 
        panel.grid.minor = element_blank(), panel.border = element_blank(), 
        panel.background = element_blank(), plot.title = element_text(size = 14, 
            face = "bold"), axis.text.x = element_text(colour = "black", 
            size = 9), axis.text.y = element_text(colour = "black", size = 9))
    
    
    df2 <- data.frame(clustcoef = cc * length(cc))
    pt2 <- ggplot(df2, aes(x = df2[["clustcoef"]])) + geom_step(stat = "ecdf", 
        col = "red", size = 1) + scale_y_continuous(labels = percent) + 
        ggtitle("Connectivity Distribution") + xlab("Number of Genes") + 
        ylab("Cumulative Proportion") + theme_bw() + theme(axis.line = element_line(size = 1, 
        colour = "black"), panel.grid.major = element_line(colour = "#d3d3d3"), 
        panel.grid.minor = element_blank(), panel.border = element_blank(), 
        panel.background = element_blank(), plot.title = element_text(size = 14, 
            face = "bold"), axis.text.x = element_text(colour = "black", 
            size = 9), axis.text.y = element_text(colour = "black", size = 9))
    
    return(ggarrange(pt1, pt2, ncol = 2))
}
