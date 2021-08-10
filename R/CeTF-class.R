###################### CeTF class

##### setClass
#' @title
#' The CeTF Class
#'
#' @description
#' The CeTF class is data storage class that stores all the
#' results from \code{\link{runAnalysis}} function.
#' 
#' @slot Data Includes the raw, tpm and norm (see \code{\link{normExp}}) data.
#' @slot DE Includes the uniquely differentially expressed genes/TFs and
#' the statistics for all genes (see \code{\link{expDiff}}).
#' @slot Input Includes input matrices for RIF (see \code{\link{RIF}}) and PCIT 
#' (see \code{\link{PCIT}}) for both conditions in \code{\link{runAnalysis}} 
#' function analysis.
#' @slot Output Includes the matrix output from RIF analysis 
#' (see \code{\link{RIF}}) and a matrix with PCIT output, and other two matrix 
#' with raw and significant adjacency (see \code{\link{PCIT}}) for both conditions 
#' inside of \code{\link{runAnalysis}} function analysis. 
#' @slot Network Network with Gene-Gene and Gene-TF interactions for both 
#' conditions (see \code{\link{PCIT}}), main TFs resulted from the complete 
#' analysis, all the TFs identified in the input data and matrix annotating all 
#' genes and TFs.
#' 
#' @return Returns an CeTF object.
#'
#' @importFrom methods new
#' @importFrom SummarizedExperiment Assays
#'
#' @name CeTF-class
#'
#' @exportClass CeTF
#'
CeTF <- setClass(Class = "CeTF", slots = c(Data = "Assays", DE = "Assays", 
    Input = "Assays", Output = "Assays", Network = "Assays"), prototype = list(Data = Assays(), 
    DE = Assays(), Input = Assays(), Output = Assays(), Network = Assays()))


###################### show
setMethod("show", "CeTF", function(object) {
    cat("class:", class(object), "\n")
    cat("Variables:", length(slotNames(object)), "\n")
    
    ## 
    nms <- slotNames(object)
    cat(sprintf("Variables Names: %s\n", paste(nms, collapse = ", ")))
})


###################### setGenerics getData
#' @title Data accessor for a CeTF class object.
#' 
#' @description The Data accessor access the raw, tpm and normalized data from 
#' \code{\link{runAnalysis}} function analysis.
#' 
#' @param x \code{\link{CeTF-class}} object
#' @param type Type of data: raw, tpm or norm (default: raw)
#' 
#' @examples
#' # load the CeTF class object resulted from runAnalysis function
#' data(CeTFdemo)
#'
#' getData(CeTFdemo)
#' 
#' @seealso \code{\link{runAnalysis}}.
#' 
#' @name getData
#' 
#' @rdname getData-methods
#' 
#' @return Returns the raw, tpm or normalized data.
#' 
#' @exportMethod getData
setGeneric(name = "getData", def = function(x, type = "raw") {
    standardGeneric("getData")
})


############# getDE
#' @title Differential Expression accessor for a CeTF class object.
#' 
#' @description The DE accessor access the differential expression resulted from 
#' \code{\link{runAnalysis}} function analysis.
#' 
#' @param x \code{\link{CeTF-class}} object
#' @param type Type of DE matrix: unique and all (default: unique)
#' 
#' @examples
#' # load the CeTF class object resulted from runAnalysis function
#' data(CeTFdemo)
#'
#' getDE(CeTFdemo)
#' 
#' @seealso \code{\link{runAnalysis}}.
#' 
#' @name getDE
#' 
#' @rdname getDE-methods
#' 
#' @return Returns the DE genes with the statistics.
#' 
#' @exportMethod getDE
setGeneric(name = "getDE", def = function(x, type = "unique") {
    standardGeneric("getDE")
})


############# InputData
#' @title Input data accessor for a CeTF class object.
#' 
#' @description The input accessor access the input matrices used for RIF and 
#' PCIT analysis to both conditions resulted from \code{\link{runAnalysis}} 
#' function analysis.
#' 
#' @param x \code{\link{CeTF-class}} object
#' @param analysis Type of analysis: rif, pcit1, pcit2. The numbers 1 and 2 
#' correspond to the respective condition (default: rif).
#' 
#' @examples
#' # load the CeTF class object resulted from runAnalysis function
#' data(CeTFdemo)
#'
#' InputData(CeTFdemo)
#' 
#' @seealso \code{\link{runAnalysis}}.
#' 
#' @name InputData
#' 
#' @rdname InputData-methods
#' 
#' @return Returns the Inputs used for RIF and PCIT.
#' 
#' @exportMethod InputData
setGeneric(name = "InputData", def = function(x, analysis = "rif") {
    standardGeneric("InputData")
})


############# OutputData
#' @title Output data accessor for a CeTF class object.
#' 
#' @description The output accessor access the output matrices and lists 
#' used for RIF and PCIT analysis to both conditions resulted from 
#' \code{\link{runAnalysis}} function analysis.
#' 
#' @param x \code{\link{CeTF-class}} object
#' @param analysis Type of analysis: rif, pcit1, pcit2. The numbers 1 and 2 
#' correspond to the respective condition (default: rif).
#' @param type Type of matrix for PCIT output: tab, adj_raw or adj_sig 
#' (default: tab).
#' 
#' @examples
#' # load the CeTF class object resulted from runAnalysis function
#' data(CeTFdemo)
#'
#' OutputData(CeTFdemo)
#' 
#' @seealso \code{\link{runAnalysis}}.
#' 
#' @name OutputData
#' 
#' @rdname OutputData-methods
#' 
#' @return Returns the Outputs used for RIF and PCIT.
#' 
#' @exportMethod OutputData
setGeneric(name = "OutputData", def = function(x, analysis = "rif", type = "tab") {
    standardGeneric("OutputData")
})


############# NetworkData
#' @title Networks data accessor for a CeTF class object.
#' 
#' @description The networks accessor access the networks, key TFs and 
#' annotations for each gene and TF resulted from PCIT analysis and
#' \code{\link{runAnalysis}} function analysis.
#' 
#' @param x \code{\link{CeTF-class}} object
#' @param type Type of data: network1, network2, keytfs, tfs or annotation.The 
#' numbers 1 and 2 correspond to the respective condition (default: network1).
#' 
#' @examples
#' # load the CeTF class object resulted from runAnalysis function
#' data(CeTFdemo)
#'
#' NetworkData(CeTFdemo)
#' 
#' @seealso \code{\link{runAnalysis}}.
#' 
#' @name NetworkData
#' 
#' @rdname NetworkData-methods
#' 
#' @return Returns the Outputs used for RIF and PCIT.
#' 
#' @exportMethod NetworkData
setGeneric(name = "NetworkData", def = function(x, type = "network1") {
    standardGeneric("NetworkData")
})





###################### setMethods getData
#' @rdname getData-methods
#' 
#' @aliases getData,CeTF-method
#' 
setMethod(f = "getData", signature = "CeTF", definition = function(x, type = "raw") {
    if (type == "raw") {
        x@Data@data$raw
    } else if (type == "tpm") {
        x@Data@data$tpm
    } else if (type == "norm") {
        x@Data@data$norm
    }
})


############# getDE
#' @rdname getDE-methods
#' 
#' @aliases getDE,CeTF-method
#' 
setMethod(f = "getDE", signature = "CeTF", definition = function(x, type = "unique") {
    if (type == "unique") {
        x@DE@data$DE_unique
    } else if (type == "all") {
        x@DE@data$DE
    }
})


############# InputData
#' @rdname InputData-methods
#' 
#' @aliases InputData,CeTF-method
#' 
setMethod(f = "InputData", signature = "CeTF", definition = function(x, 
    analysis = "rif") {
    if (analysis == "rif") {
        x@Input@data$RIF_input
    } else if (analysis == "pcit1") {
        x@Input@data$PCIT_input_cond1
    } else if (analysis == "pcit2") {
        x@Input@data$PCIT_input_cond2
    }
})


############# OutputData
#' @rdname OutputData-methods
#' 
#' @aliases OutputData,CeTF-method
#' 
setMethod(f = "OutputData", signature = "CeTF", definition = function(x, 
    analysis = "rif", type = "tab") {
    if (analysis == "rif") {
        x@Output@data$RIF_out
    } else if (analysis == "pcit1" & type == "tab") {
        x@Output@data$PCIT_out_cond1$tab
    } else if (analysis == "pcit1" & type == "adj_raw") {
        x@Output@data$PCIT_out_cond1$adj_raw
    } else if (analysis == "pcit1" & type == "adj_sig") {
        x@Output@data$PCIT_out_cond1$adj_sig
    } else if (analysis == "pcit2" & type == "tab") {
        x@Output@data$PCIT_out_cond2$tab
    } else if (analysis == "pcit2" & type == "adj_raw") {
        x@Output@data$PCIT_out_cond2$adj_raw
    } else if (analysis == "pcit2" & type == "adj_sig") {
        x@Output@data$PCIT_out_cond2$adj_sig
    }
})


############# NetworkData
#' @rdname NetworkData-methods
#' 
#' @aliases NetworkData,CeTF-method
#' 
setMethod(f = "NetworkData", signature = "CeTF", definition = function(x, 
    type = "network1") {
    if (type == "network1") {
        x@Network@data$net_cond1
    } else if (type == "network2") {
        x@Network@data$net_cond2
    } else if (type == "keytfs") {
        x@Network@data$keyTF
    } else if (type == "tfs") {
        x@Network@data$tfs
    } else if (type == "annotation") {
        x@Network@data$anno
    }
})


### Changing Validity Assays to accept objects of different dims
#' @export
#' @importFrom  S4Vectors setValidity2
.valid.Assays <- function(x) {
    assays <- as(x, "SimpleList", strict = FALSE)
    if (!is(assays, "SimpleList")) 
        return("'assays' must be a SimpleList object")
    if (length(assays) == 0L) 
        return(NULL)
    
    NULL
}
setValidity2("Assays", .valid.Assays)
