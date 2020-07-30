#' @title Circos plot for the Transcription Factors/genes targets.
#' 
#' @description Generate an plot for Transcription Targets (TFs) or any gene 
#' targets. This plot consists of sorting all the chromosomes of any specie
#' based in GTF annotation file and showing how the selected TF(s)/gene(s) 
#' targets are distributed. If a target is connected to the same chromosome as 
#' the selected one so the connection is defined as \emph{cis}, otherwise it is
#' a \emph{trans} connection.
#' 
#' @param object CeTF class object resulted from \code{\link{runAnalysis}} 
#' function.
#' @param file GTF file or path.
#' @param nomenclature Gene nomenclature: \emph{SYMBOL} or \emph{ENSEMBL.}
#' @param selection Specify a single or up to 4 TF/gene to be visualized for.
#' @param cond The options are \emph{condition1} or \emph{condition2} based on
#' the conditions previously defined in \code{\link{runAnalysis}} function.
#' 
#' @return Returns an plot with a specific(s) TF/gene and its targets in order
#' to visualize the chromosome location of each one.
#'
#' @details The black links are between different chromosomes while the red 
#' links are between the same chromosome.
#'
#' @importFrom GenomicTools gtfToBed
#' @importFrom GenomicTools.fileHandler importGTF
#' @importFrom circlize circos.par circos.initialize circos.track circos.rect 
#' rand_color circos.text circos.link highlight.chromosome circos.clear CELL_META 
#' @importFrom graphics title layout
#'
#' @examples
#' \dontrun{
#' CircosTargets(object = out, 
#' file = '/path/to/gtf/specie.gtf', 
#' nomenclature = 'SYMBOL', 
#' selection = 'TCF4', 
#' cond = 'condition1')
#' }
#'
#' @export
CircosTargets <- function(object, file, nomenclature, selection, cond) {
    if (!is(object, "CeTF")) {
        stop("the input must be a CeTF class object")
    }
    if (missing(object)) {
        stop("No \"object\" parameter provided")
    }
    if (missing(file)) {
        stop("No \"file\" parameter provided")
    }
    
    if (cond == "condition1") {
        net <- NetworkData(object, "network1")
    } else if (cond == "condition2") {
        net <- NetworkData(object, "network2")
    } else {
        stop("Input a valid condition: condition1 or condition2")
    }
    
    genes <- unique(c(as.character(net[which(net[["gene1"]] %in% selection), 
        2]), as.character(net[which(net[["gene2"]] %in% selection), 1])))
    
    if (length(genes) == 0) {
        stop("This/These gene(s) are not on the network.")
    }
    
    GTF <- importGTF(file, verbose = FALSE)
    annotBed <- gtfToBed(GTF)
    
    if (nomenclature == "SYMBOL") {
        bedOut <- data.frame(Chr = GTF[, 1], Start = GTF[, 4], End = GTF[, 
            5], Gene = GTF[["gene_name"]], stringsAsFactors = FALSE)
    } else if (nomenclature == "ENSEMBL") {
        bedOut <- data.frame(Chr = GTF[, 1], Start = GTF[, 4], End = GTF[, 
            5], Gene = GTF[["gene_id"]], stringsAsFactors = FALSE)
    } else {
        stop("nomenclature should be SYMBOL or ENSEMBL")
    }
    colnames(bedOut) <- c("chr", "start", "end", "gene")
    
    annotBed <- subset(bedOut, bedOut[["chr"]] %in% c(seq_len(50), "X", 
        "Y"))
    
    convert <- unique(annotBed[["chr"]][order(annotBed[["chr"]])])
    numbers <- gsub("[^0-9.]", "", convert)
    numbers <- as.numeric(numbers[numbers != ""])
    letters <- sub("^([[:alpha:]]*).*", "\\1", convert)
    letters <- letters[letters != ""]
    annotBed[["chr"]] <- factor(annotBed[["chr"]], levels = c(sort(numbers), 
        sort(letters)))
    
    chr <- unique(annotBed[["chr"]])
    ln <- NULL
    for (i in seq_along(chr)) {
        tmp1 <- annotBed[which(annotBed[["chr"]] == chr[i]), ]
        ln <- rbind(ln, c(as.character(chr[i]), max(tmp1[["end"]]) - min(tmp1[["start"]])))
    }
    
    if (length(selection) == 1) {
        layout(matrix(1, 1, 1))
        circos.par(cell.padding = c(0.02, 0, 0.02, 0))
        circos.initialize(annotBed[["chr"]], x = annotBed[["start"]], xlim = c(0, 
            1), sector.width = as.numeric(ln[, 2]))
        circos.track(annotBed[["chr"]], ylim = c(0, 1), panel.fun = function(x, 
            y) {
            chr = CELL_META$sector.index
            xlim = CELL_META$xlim
            ylim = CELL_META$ylim
            circos.rect(xlim[1], 0, xlim[2], 1, col = rand_color(1))
            circos.text(mean(xlim), mean(ylim), chr, cex = 1, col = "white", 
                facing = "inside", niceFacing = TRUE)
        }, track.height = 0.15, bg.border = NA)
        
        bed2 <- annotBed[annotBed[["gene"]] %in% genes, ]
        bed2 <- bed2[!duplicated(bed2[["gene"]]), ]
        bed2[["value1"]] <- bed2[["end"]] - bed2[["start"]]
        bed2[["value1"]] <- (bed2[["value1"]] - min(bed2[["value1"]]))/(max(bed2[["value1"]]) - 
            min(bed2[["value1"]]))
        bed2 <- bed2[, c(1, 2, 3, 5)]
        
        bed1 <- annotBed[annotBed[["gene"]] == selection, ]
        bed1 <- bed1[!duplicated(bed1[["gene"]]), ]
        bed1[["value1"]] <- (((bed1[["end"]] - bed1[["start"]]) * 100)/bed1[["end"]])
        bed1 <- do.call("rbind", replicate(nrow(bed2), bed1[1, ], simplify = FALSE))
        bed1 <- bed1[, c(1, 2, 3, 5)]
        
        for (i in seq_along(bed1[["chr"]])) {
            if (bed1[["chr"]][i] == bed2[["chr"]][i]) {
                circos.link(bed1[["chr"]][i], bed1[["value1"]][i], bed2[["chr"]][i], 
                  bed2[["value1"]][i], border = 1, col = "red")
            } else {
                circos.link(bed1[["chr"]][i], bed1[["value1"]][i], bed2[["chr"]][i], 
                  bed2[["value1"]][i], border = 1, col = "#00000080")
            }
        }
        
        highlight.chromosome(unique(bed1[["chr"]]))
        title(paste("Circos Plot for", selection))
        circos.clear()
    } else if (length(selection) > 1) {
        layout(matrix(seq_len(4), 2, 2))
        for (i in seq_along(selection)) {
            circos.par(cell.padding = c(0.02, 0, 0.02, 0))
            circos.initialize(annotBed[["chr"]], x = annotBed[["start"]], 
                xlim = c(0, 1), sector.width = as.numeric(ln[, 2]))
            
            circos.track(annotBed[["chr"]], ylim = c(0, 1), panel.fun = function(x, 
                y) {
                chr = CELL_META$sector.index
                xlim = CELL_META$xlim
                ylim = CELL_META$ylim
                circos.rect(xlim[1], 0, xlim[2], 1, col = rand_color(1))
                circos.text(mean(xlim), mean(ylim), chr, cex = 0.7, col = "white", 
                  facing = "inside", niceFacing = TRUE)
            }, track.height = 0.15, bg.border = NA)
            
            genes <- c(as.character(subset(net, net[["gene1"]] == selection[i])[, 
                2]), as.character(subset(net, net[["gene2"]] == selection[i])[, 
                1]))
            
            bed2 <- annotBed[annotBed[["gene"]] %in% genes, ]
            bed2 <- bed2[!duplicated(bed2[["gene"]]), ]
            bed2[["value1"]] <- bed2[["end"]] - bed2[["start"]]
            bed2[["value1"]] <- (bed2[["value1"]] - min(bed2[["value1"]]))/(max(bed2[["value1"]]) - 
                min(bed2[["value1"]]))
            bed2 <- bed2[, c(1, 2, 3, 5)]
            
            bed1 <- annotBed[annotBed[["gene"]] == selection[i], ]
            bed1 <- bed1[!duplicated(bed1[["gene"]]), ]
            bed1[["value1"]] <- (((bed1[["end"]] - bed1[["start"]]) * 100)/bed1[["end"]])
            bed1 <- do.call("rbind", replicate(nrow(bed2), bed1[1, ], simplify = FALSE))
            bed1 <- bed1[, c(1, 2, 3, 5)]
            
            for (j in seq_along(bed1[["chr"]])) {
                if (bed1[["chr"]][j] == bed2[["chr"]][j]) {
                  circos.link(bed1[["chr"]][j], bed1[["value1"]][j], bed2[["chr"]][j], 
                    bed2[["value1"]][j], border = 1, col = "red")
                } else {
                  circos.link(bed1[["chr"]][j], bed1[["value1"]][j], bed2[["chr"]][j], 
                    bed2[["value1"]][j], border = 1, col = "#00000080")
                }
            }
            highlight.chromosome(unique(bed1[["chr"]]))
            title(selection[i])
            circos.clear()
        }
    }
}
