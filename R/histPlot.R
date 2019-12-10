#' @title Histogram of connectivity distribution
#'
#' @description Generate the histogram for adjacency matrix to show the clustering coefficient distribution.
#'
#' @param mat Adjacency matrix resulting from PCIT analysis in which has some zero values.
#'
#' @return The histogram of connectivity distribution.
#'
#' @importFrom ggplot2 ggplot geom_histogram xlab ylab element_line element_blank
#' @importFrom graphics hist
#' @importFrom gridExtra grid.arrange
#'
#' @examples
#' data('simNorm')
#' results <- PCIT(simNorm)
#' histPlot(results[[3]])
#'
#'
#'
#' @export
histPlot <- function(mat) {
    if (!is.data.frame(mat) & !is.matrix(mat)) {
        stop("input must be a dataframe or a matrix")
    }
    
    cc <- clustCoef(mat)
    
    df1 <- data.frame(clustcoef = cc * length(cc))
    pt1 <- ggplot(df1, aes(clustcoef)) + geom_histogram(breaks = seq(0, 100, by = 10), col = "black", 
        fill = "#1F3552") + ggtitle("Connectivity Distribution") + xlab("Number of Connections") + 
        ylab("Number of Genes") + scale_y_continuous(expand = c(0, 0)) + theme_bw() + theme(axis.line = element_line(size = 1, 
        colour = "black"), panel.grid.major = element_line(colour = "#d3d3d3"), panel.grid.minor = element_blank(), 
        panel.border = element_blank(), panel.background = element_blank(), plot.title = element_text(size = 14, 
            face = "bold"), axis.text.x = element_text(colour = "black", size = 9), axis.text.y = element_text(colour = "black", 
            size = 9))
    
    df2 <- data.frame(clustcoef = cc)
    pt2 <- ggplot(df2, aes(clustcoef)) + geom_histogram(breaks = seq(0, 0.6, by = 0.05), col = "black", 
        fill = "#1F3552") + ggtitle("Connectivity Distribution") + xlab("Proportion of Connections") + 
        ylab("Number of Genes") + scale_y_continuous(expand = c(0, 0)) + theme_bw() + theme(axis.line = element_line(size = 1, 
        colour = "black"), panel.grid.major = element_line(colour = "#d3d3d3"), panel.grid.minor = element_blank(), 
        panel.border = element_blank(), panel.background = element_blank(), plot.title = element_text(size = 14, 
            face = "bold"), axis.text.x = element_text(colour = "black", size = 9), axis.text.y = element_text(colour = "black", 
            size = 9))
    
    return(grid.arrange(grobs = list(pt1, pt2), ncol = 2))
}
