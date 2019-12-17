#' @title Density Plot of raw correlation coefficients
#'
#' @description Generate the density plot for an adjacency matrix.
#'
#' @param mat An raw adjacency matrix.
#'
#' @return A density plot for raw correlation coefficients.
#'
#' @importFrom ggplot2 ggplot aes geom_density scale_x_continuous scale_y_continuous ggtitle theme_bw theme geom_vline element_text element_line element_blank
#'
#' @examples
#' # loading a simulated normalized data
#' data(simNorm)
#'
#' # getting the PCIT results
#' results <- PCIT(simNorm)
#'
#' # using the PCIT results, more specifically the raw adjacency
#' matrix to get the density plot
#' densityPlot(results$adj_raw)
#'
#' @export
densityPlot <- function(mat) {
    if (!is.matrix(mat)) {
        stop("mat must be a dataframe or a matrix")
    }

    df <- data.frame(corr = mat[upper.tri(mat)])
    pt <- ggplot(df, aes(x = corr)) + geom_density(colour = "black",
        fill = "#a0b8d6", size = 1) + scale_x_continuous(name = "Correlation Coefficient",
        breaks = seq(-1, 1, 0.2), limits = c(-1.15,
            1.15)) + scale_y_continuous(name = "Density",
        expand = c(0, 0), limits = c(0, 1.1)) + ggtitle("Density Plot of Raw Correlation Coefficients") +
        theme_bw() + theme(axis.line = element_line(size = 1,
        colour = "black"), panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), plot.title = element_text(size = 14,
            face = "bold"), axis.text.x = element_text(colour = "black",
            size = 9), axis.text.y = element_text(colour = "black",
            size = 9)) + geom_vline(xintercept = mean(df$corr),
        size = 1, colour = "#FF3721", linetype = "dashed")
    return(pt)
}
