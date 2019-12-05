#' @title density plot
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
density.sig <- function(mat1, mat2, type = NULL) {
  if(is.null(type)){stop("a type of plot is needed")}

  if (type == 1) {
    df1 <- data.frame(corr = mat1[upper.tri(mat1) & !is.na(mat1)],
                      sig = "All")
    df2 <- data.frame(corr = mat2[intersect(which(as.single(mat2) != 0), which(upper.tri(mat1)))],
                      sig = "PCIT Significant")
    df <- rbind(df1, df2)

    pt <- ggplot(df,aes(corr))+
      geom_histogram(data=subset(df,sig=='All'),aes(fill=sig), breaks=seq(-1, 1, by = 0.05), col="black")+
      geom_histogram(data=subset(df,sig=='PCIT Significant'),aes(fill=sig), breaks=seq(-1, 1, by = 0.05), col="black")+
      scale_fill_manual(name="", values=c("gray","black"),labels=c("All","PCIT Significant")) +
      ggtitle("Density Distribution of Correlation Coefficients") +
      xlab("Correlation Coefficient") + ylab("Frequency") +
      scale_y_continuous(expand = c(0, 0), limits=c(0,mean(hist(df$corr, plot = F)$counts)+50)) +
      theme_bw() +
      theme(axis.line = element_line(size=1, colour = "black"),
            panel.grid.major = element_line(colour = "#d3d3d3"),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(), panel.background = element_blank(),
            plot.title = element_text(size = 14, face = "bold"),
            axis.text.x=element_text(colour="black", size = 9),
            axis.text.y=element_text(colour="black", size = 9),
            legend.position="top")
  } else if (type == 2) {
    df1 <- data.frame(corr = mat1[upper.tri(mat1) & !is.na(mat1)],
                      sig = "All")
    df2 <- data.frame(corr = mat1[intersect(which(as.single(mat2) != 0), which(upper.tri(mat2)))],
                      sig = "PCIT Significant")
    df3 <- data.frame(corr = mat1[intersect(which(abs(mat1) > 0.5), which(upper.tri(mat1)))],
                      sig = "abs. cor. > 0.5")
    df <- rbind(df1, df2, df3)

    pt <- ggplot(df,aes(corr))+
      geom_histogram(data=subset(df,sig=='All'),aes(fill=sig), breaks=seq(-1, 1, by = 0.05), col="black") +
      geom_histogram(data=subset(df,sig=='abs. cor. > 0.5'),aes(fill=sig), breaks=seq(-1, 1, by = 0.05), col="black") +
      geom_histogram(data=subset(df,sig=='PCIT Significant'),aes(fill=sig), breaks=seq(-1, 1, by = 0.05), col="black") +
      scale_fill_manual(name="", values=c("red","gray","black"),labels=c("abs. cor. > 0.5", "All", "PCIT Significant")) +
      ggtitle("Density Distribution of Correlation Coefficients") +
      xlab("Correlation Coefficient") + ylab("Frequency") +
      scale_y_continuous(expand = c(0, 0), limits=c(0,mean(hist(df$corr, plot = F)$counts)+50)) +
      theme_bw() +
      theme(axis.line = element_line(size=1, colour = "black"),
            panel.grid.major = element_line(colour = "#d3d3d3"),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(), panel.background = element_blank(),
            plot.title = element_text(size = 14, face = "bold"),
            axis.text.x=element_text(colour="black", size = 9),
            axis.text.y=element_text(colour="black", size = 9),
            legend.position="top")
  } else if (type == 3) {
    df1 <- data.frame(corr = mat1[upper.tri(mat1) & !is.na(mat1)],
                      sig = "All")
    df2 <- data.frame(corr = mat1[intersect(which(as.single(mat2) != 0), which(upper.tri(mat2)))],
                      sig = "PCIT Significant")
    df3 <- data.frame(corr = mat1[intersect(which(abs(mat1) > 0.5), which(upper.tri(mat2)))],
                      sig = "abs. cor. > 0.5")
    df <- rbind(df1, df2, df3)

    pt <- ggplot(df,aes(corr))+
      geom_histogram(data=subset(df,sig=='All'),aes(fill=sig), breaks=seq(-1, 1, by = 0.05), col="black", alpha = 0.8) +
      geom_histogram(data=subset(df,sig=='PCIT Significant'),aes(fill=sig), breaks=seq(-1, 1, by = 0.05), col="black", alpha = 0.9) +
      geom_histogram(data=subset(df,sig=='abs. cor. > 0.5'),aes(fill=sig), breaks=seq(-1, 1, by = 0.05), col="black", alpha = 0.8) +
      scale_fill_manual(name="", values=c("gray10","gray","darkred"),labels=c("abs. cor. > 0.5", "All", "PCIT Significant")) +
      ggtitle("Density Distribution of Correlation Coefficients") +
      xlab("Correlation Coefficient") + ylab("Frequency") +
      scale_y_continuous(expand = c(0, 0), limits=c(0,mean(hist(df$corr, plot = F)$counts)+50)) +
      theme_bw() +
      theme(axis.line = element_line(size=1, colour = "black"),
            panel.grid.major = element_line(colour = "#d3d3d3"),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(), panel.background = element_blank(),
            plot.title = element_text(size = 14, face = "bold"),
            axis.text.x=element_text(colour="black", size = 9),
            axis.text.y=element_text(colour="black", size = 9),
            legend.position="top")
  }
  return(pt)
}
