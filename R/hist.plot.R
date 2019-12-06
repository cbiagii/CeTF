#' @title hist plot
#'
#' @description teste.
#'
#' @param mat teste
#'
#' @return teste.
#'
#' @importFrom ggplot2 ggplot
#'
#' @examples
#' teste
#'
#' @export
hist.plot <- function(mat, type = c("number", "prop")) {
  if (!is.data.frame(mat) & !is.matrix(mat)) {
    stop("input must be a dataframe or a matrix")
  }

  type <- match.arg(type, choices=c("number", "prop"))

  cc <- clustCoef(mat)

  if (type == "number") {
    df <- data.frame(clustcoef = cc*length(cc))
    pt <- ggplot(df, aes(clustcoef)) +
      geom_histogram(breaks=seq(0, 100, by = 10), col="black", fill="#1F3552") +
      ggtitle("Connectivity Distribution") +
      xlab("Number of Connections") + ylab("Number of Genes") +
      scale_y_continuous(expand = c(0, 0), limits=c(0,max(hist(df$clustcoef, plot = F)$counts)+2)) +
      theme_bw() +
      theme(axis.line = element_line(size=1, colour = "black"),
            panel.grid.major = element_line(colour = "#d3d3d3"),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(), panel.background = element_blank(),
            plot.title = element_text(size = 14, face = "bold"),
            axis.text.x=element_text(colour="black", size = 9),
            axis.text.y=element_text(colour="black", size = 9))
  } else if (type == "prop") {
    df <- data.frame(clustcoef = cc)

    pt <- ggplot(df, aes(clustcoef)) +
      geom_histogram(breaks=seq(0, 0.6, by = 0.05), col="black", fill="#1F3552") +
      ggtitle("Connectivity Distribution") +
      xlab("Proportion of Connections") + ylab("Number of Genes") +
      scale_y_continuous(expand = c(0, 0), limits=c(0,max(hist(df$clustcoef, plot = F)$counts)+2)) +
      theme_bw() +
      theme(axis.line = element_line(size=1, colour = "black"),
            panel.grid.major = element_line(colour = "#d3d3d3"),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(), panel.background = element_blank(),
            plot.title = element_text(size = 14, face = "bold"),
            axis.text.x=element_text(colour="black", size = 9),
            axis.text.y=element_text(colour="black", size = 9))
  }
  return(pt)
}
