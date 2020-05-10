#' @title Relationship plots between RIF1, RIF2 and DE genes
#'
#' @description Generate plots for the relationship between the RIF output 
#' analysis (RIF1 and RIF2) and for differentially expressed genes (DE).
#'
#' @param object CeTF object resulted from \code{\link{runAnalysis}} function.
#' @param color Color of points (default: darkblue)
#' @param type Type of plot. The available options are: RIF or 
#' DE (default: RIF)
#'
#' @return Returns a relationship plot between RIF1 and RIF2 or a plot with the 
#' relationship between RIF1 or RIF2 with DE genes.
#'
#' @importFrom ggpubr ggarrange
#' @importFrom ggplot2 ggplot geom_point xlim ylim xlab ylab theme_bw
#'
#' @details 
#' This function can only be used after using the \code{\link{runAnalysis}} 
#' function as it uses the CeTF class object as input.
#' 
#' @examples
#' # load the CeTF class object resulted from runAnalysis function
#' data(CeTFdemo)
#'
#' # performing RIFPlot for RIF
#' RIFPlot(object = CeTFdemo, 
#'         color  = 'darkblue', 
#'         type   = 'RIF')
#' 
#'         
#' # performing RIFPlot for DE
#' RIFPlot(object = CeTFdemo, 
#'         color  = 'darkblue', 
#'         type   = 'DE')
#'
#' @export
RIFPlot <- function(object, color = "darkblue", type = "RIF") {
    if (!is(object, "CeTF")) {
        stop("the input must be a CeTF class object")
    }
    
    if (type == "RIF") {
        tmp <- OutputData(object, analysis = "rif")
        
        pt <- ggplot(data = tmp, aes(x = .data[["RIF1"]], y = .data[["RIF2"]])) + 
            geom_point(color = color, size = 1.5) + xlim(-6, 6) + ylim(-6, 
            6) + xlab("RIF1") + ylab("RIF2") + theme_bw()
        
        return(pt)
    } else if (type == "DE") {
        tmp1 <- OutputData(object, analysis = "rif")
        
        tmp2 <- getDE(object, "all")
        tmp2 <- tmp2[rownames(tmp2) %in% tmp1$TF, ]
        
        tb <- data.frame(RIF1 = tmp1$RIF1, RIF2 = tmp1$RIF2, DE = tmp2$diff)
        
        pt1 <- ggplot(data = tb, aes(x = .data[["RIF1"]], y = .data[["DE"]])) + 
            geom_point(color = color, size = 1.5) + xlab("RIF1") + ylab("Expression Difference") + 
            theme_bw()
        
        pt2 <- ggplot(data = tb, aes(x = .data[["RIF2"]], y = .data[["DE"]])) + 
            geom_point(color = color, size = 1.5) + xlab("RIF2") + ylab("Expression Difference") + 
            theme_bw()
        
        pt3 <- ggarrange(pt1, pt2, nrow = 2, ncol = 1)
        
        return(pt3)
    }
}
