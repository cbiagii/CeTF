#' @title Converts GTF to BED
#'
#' @description Converts a GTF to BED format.
#'
#' @param gtf A GTF as data.frame.
#'
#' @return Returns a data.frame in BED format.
#'
#' @export
gtfToBed <- function(gtf){
  correction <- 0
  if(sum(class(gtf) == "gtf") == 0) stop("This function requires a gtf-object as import (as given by the importGTF() function)")
  bedOut <- data.frame(Chr = gtf[,1],
                       Start = gtf[,4] - 1,
                       End = gtf[,5] - 1,
                       Gene = gtf$gene_id,
                       stringsAsFactors = FALSE)
  
  colnames(bedOut) <- c("Chr", "Start", "End", "Gene")
  class(bedOut) <- append(class(bedOut), "bed")
  bedOut
}