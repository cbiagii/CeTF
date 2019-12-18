#' @title Density distribution of correlation coefficients and significant PCIT values
#'
#' @description Generate the density plot for adjacency matrices. This
#' function uses the raw adjacency matrix and significant adjacency
#' matrix resulted from \code{\link{PCIT}} function.
#'
#' @param mat1 Raw adjacency matrix.
#' @param mat2 Significant adjacency matrix.
#' @param threshold Threshold of correlation module to plot (default: 0.5).
#'
#' @return a density plot of raw correlation with significant PCIT values.
#'
#' @importFrom ggplot2 geom_histogram scale_fill_manual xlab ylab element_line element_blank
#' @importFrom ggpubr ggarrange
#'
#' @examples
#' # loading a simulated normalized data
#' data('simNorm')
#'
#' # getting the PCIT results
#' results <- PCIT(simNorm)
#'
#' # using the PCIT results to get density distribution of correlation coefficients
#' densitySig(mat1 = results$adj_raw,
#'            mat2 = results$adj_sig,
#'            threshold = 0.5)
#'
#' @export
densitySig <- function(mat1, mat2, threshold = 0.5) {
<<<<<<< HEAD
    if(!is.data.frame(mat1) & !is.matrix(mat1)){stop("mat1 must be a dataframe or a matrix")}
    if(!is.data.frame(mat2) & !is.matrix(mat2)){stop("mat2 must be a dataframe or a matrix")}

    df1 <- data.frame(corr = mat1[upper.tri(mat1) &
=======
    if (!is.data.frame(mat1) & !is.matrix(mat1)) {
        stop("mat1 must be a dataframe or a matrix")
    }
    if (!is.data.frame(mat2) & !is.matrix(mat2)) {
        stop("mat2 must be a dataframe or a matrix")
    }
    
    df1 <- data.frame(corr = mat1[upper.tri(mat1) & 
>>>>>>> 43614a53fc5fd047595c36314fe49c8a0a0915a2
        !is.na(mat1)], sig = "All")
    df2 <- data.frame(corr = mat2[intersect(which(as.single(mat2) != 
        0), which(upper.tri(mat1)))], sig = "PCIT Significant")
    df3 <- data.frame(corr = mat1[intersect(which(abs(mat1) > 
        threshold), which(upper.tri(mat1)))], sig = paste0("abs. cor. > ", 
        threshold))
    
    df <- rbind(df1, df2)
    pt1 <- ggplot(df, aes(corr)) + geom_histogram(data = subset(df, 
        sig == "All"), aes(fill = sig), breaks = seq(-1, 
        1, by = 0.05), col = "black") + geom_histogram(data = subset(df, 
        sig == "PCIT Significant"), aes(fill = sig), 
        breaks = seq(-1, 1, by = 0.05), col = "black") + 
        scale_fill_manual(name = "", values = c("gray", 
            "black"), labels = c("All", "PCIT Significant")) + 
        ggtitle("Density Distribution of Correlation Coefficients") + 
        xlab("Correlation Coefficient") + ylab("Frequency") + 
        scale_y_continuous(expand = c(0, 0)) + theme_bw() + 
        theme(axis.line = element_line(size = 1, colour = "black"), 
            panel.grid.major = element_line(colour = "#d3d3d3"), 
            panel.grid.minor = element_blank(), panel.border = element_blank(), 
            panel.background = element_blank(), plot.title = element_text(size = 14, 
                face = "bold"), axis.text.x = element_text(colour = "black", 
                size = 9), axis.text.y = element_text(colour = "black", 
                size = 9), legend.position = "top")
    
    df <- rbind(df1, df2, df3)
    pt2 <- ggplot(df, aes(corr)) + geom_histogram(data = subset(df, 
        sig == "All"), aes(fill = sig), breaks = seq(-1, 
        1, by = 0.05), col = "black") + geom_histogram(data = subset(df, 
        sig == paste0("abs. cor. > ", threshold)), 
        aes(fill = sig), breaks = seq(-1, 1, by = 0.05), 
        col = "black") + geom_histogram(data = subset(df, 
        sig == "PCIT Significant"), aes(fill = sig), 
        breaks = seq(-1, 1, by = 0.05), col = "black") + 
        scale_fill_manual(name = "", values = c("red", 
            "gray", "black"), labels = c(paste0("abs. cor. > ", 
            threshold), "All", "PCIT Significant")) + 
        xlab("Correlation Coefficient") + ylab("Frequency") + 
        scale_y_continuous(expand = c(0, 0)) + theme_bw() + 
        theme(axis.line = element_line(size = 1, colour = "black"), 
            panel.grid.major = element_line(colour = "#d3d3d3"), 
            panel.grid.minor = element_blank(), panel.border = element_blank(), 
            panel.background = element_blank(), plot.title = element_text(size = 14, 
                face = "bold"), axis.text.x = element_text(colour = "black", 
                size = 9), axis.text.y = element_text(colour = "black", 
                size = 9), legend.position = "top")
    
    df <- rbind(df1, df2, df3)
    pt3 <- ggplot(df, aes(corr)) + geom_histogram(data = subset(df, 
        sig == "All"), aes(fill = sig), breaks = seq(-1, 
        1, by = 0.05), col = "black", alpha = 0.8) + 
        geom_histogram(data = subset(df, sig == "PCIT Significant"), 
            aes(fill = sig), breaks = seq(-1, 1, by = 0.05), 
            col = "black", alpha = 0.9) + geom_histogram(data = subset(df, 
        sig == paste0("abs. cor. > ", threshold)), 
        aes(fill = sig), breaks = seq(-1, 1, by = 0.05), 
        col = "black", alpha = 0.8) + scale_fill_manual(name = "", 
        values = c("gray10", "gray", "darkred"), labels = c(paste0("abs. cor. > ", 
            threshold), "All", "PCIT Significant")) + 
        xlab("Correlation Coefficient") + ylab("Frequency") + 
        scale_y_continuous(expand = c(0, 0)) + theme_bw() + 
        theme(axis.line = element_line(size = 1, colour = "black"), 
            panel.grid.major = element_line(colour = "#d3d3d3"), 
            panel.grid.minor = element_blank(), panel.border = element_blank(), 
            panel.background = element_blank(), plot.title = element_text(size = 14, 
                face = "bold"), axis.text.x = element_text(colour = "black", 
                size = 9), axis.text.y = element_text(colour = "black", 
                size = 9), legend.position = "top")
    
    return(ggarrange(pt1, ggarrange(pt2, pt3, ncol = 2), 
        nrow = 2))
}
