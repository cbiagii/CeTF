#' The pcitrif Class
#'
#' The pcitrif class is data storage class that stores all the results of complete RIF and PCIT analysis.
#'
#' @slot step1 List of objects related to step 1 (TPM Filter) of analysis.
#' @slot step2 List of objects related to step 2 (Differential Expression) of analysis.
#' @slot step3 List of objects related to step 3 (Regulatory Impact Factors analysis) of analysis.
#' @slot step4 List of objects related to step 4 (Partial Correlation and Information Theory analysis) of analysis.
#'
#' @importFrom methods new
#' @name pcitrif-class
#' @exportClass pcitrif
#'
pcitrif <- setClass(Class = "pcitrif", slots = list(step1 = "list", step2 = "list", step3 = "list",
    step4 = "list"))



### setGenerics
#' Method rawData.
#' @name pcitrif-class
#' @rdname pcitrif-class
#' @exportMethod rawData
setGeneric("rawData", function(x) standardGeneric("rawData"))

#' Method tpmData.
#' @name pcitrif-class
#' @rdname pcitrif-class
#' @exportMethod tpmData
setGeneric("tpmData", function(x) standardGeneric("tpmData"))

#' Method getSelectedGenes.
#' @name pcitrif-class
#' @rdname pcitrif-class
#' @exportMethod getSelectedGenes
setGeneric("getSelectedGenes", function(x) standardGeneric("getSelectedGenes"))

#' Method normData.
#' @name pcitrif-class
#' @rdname pcitrif-class
#' @exportMethod normData
setGeneric("normData", function(x) standardGeneric("normData"))

#' Method getDE.
#' @name pcitrif-class
#' @rdname pcitrif-class
#' @exportMethod getDE
setGeneric("getDE", function(x) standardGeneric("getDE"))

#' Method getTF.
#' @name pcitrif-class
#' @rdname pcitrif-class
#' @exportMethod getTF
setGeneric("getTF", function(x) standardGeneric("getTF"))

#' Method getRIFinput.
#' @name pcitrif-class
#' @rdname pcitrif-class
#' @exportMethod getRIFinput
setGeneric("getRIFinput", function(x) standardGeneric("getRIFinput"))

#' Method getRIFoutput.
#' @name pcitrif-class
#' @rdname pcitrif-class
#' @exportMethod getRIFoutput
setGeneric("getRIFoutput", function(x) standardGeneric("getRIFoutput"))

#' Method getPCITinput_cond1.
#' @name pcitrif-class
#' @rdname pcitrif-class
#' @exportMethod getPCITinput_cond1
setGeneric("getPCITinput_cond1", function(x) standardGeneric("getPCITinput_cond1"))

#' Method getPCITinput_cond2.
#' @name pcitrif-class
#' @rdname pcitrif-class
#' @exportMethod getPCITinput_cond2
setGeneric("getPCITinput_cond2", function(x) standardGeneric("getPCITinput_cond2"))

#' Method getPCIToutput_cond1.
#' @name pcitrif-class
#' @rdname pcitrif-class
#' @exportMethod getPCIToutput_cond1
setGeneric("getPCIToutput_cond1", function(x) standardGeneric("getPCIToutput_cond1"))

#' Method getPCIToutput_cond2.
#' @name pcitrif-class
#' @rdname pcitrif-class
#' @exportMethod getPCIToutput_cond2
setGeneric("getPCIToutput_cond2", function(x) standardGeneric("getPCIToutput_cond2"))

#' Method getNetworkCond1.
#' @name pcitrif-class
#' @rdname pcitrif-class
#' @exportMethod getNetworkCond1
setGeneric("getNetworkCond1", function(x) standardGeneric("getNetworkCond1"))

#' Method getNetworkCond2.
#' @name pcitrif-class
#' @rdname pcitrif-class
#' @exportMethod getNetworkCond2
setGeneric("getNetworkCond2", function(x) standardGeneric("getNetworkCond2"))

#' Method getRIFkeyTF.
#' @name pcitrif-class
#' @rdname pcitrif-class
#' @exportMethod getRIFkeyTF
setGeneric("getRIFkeyTF", function(x) standardGeneric("getRIFkeyTF"))

#' Method getAnno.
#' @name pcitrif-class
#' @rdname pcitrif-class
#' @exportMethod getAnno
setGeneric("getAnno", function(x) standardGeneric("getAnno"))



# setMethods
#' @rdname pcitrif-class
#' @aliases rawData,pcitrif-method
setMethod("rawData", "pcitrif", function(x) x@step1$raw)

#' @rdname pcitrif-class
#' @aliases tpmDat,pcitrif-method
setMethod("tpmData", "pcitrif", function(x) x@step1$tpm)

#' @rdname pcitrif-class
#' @aliases getSelectedGenes,pcitrif-method
setMethod("getSelectedGenes", "pcitrif", function(x) x@step1$selected_genes)

#' @rdname pcitrif-class
#' @aliases normData,pcitrif-method
setMethod("normData", "pcitrif", function(x) x@step1$norm)

#' @rdname pcitrif-class
#' @aliases getDE,pcitrif-method
setMethod("getDE", "pcitrif", function(x) x@step2$de_genes)

#' @rdname pcitrif-class
#' @aliases getTF,pcitrif-method
setMethod("getTF", "pcitrif", function(x) x@step2$tf)

#' @rdname pcitrif-class
#' @aliases getRIFinput,pcitrif-method
setMethod("getRIFinput", "pcitrif", function(x) x@step3$input)

#' @rdname pcitrif-class
#' @aliases getRIFoutput,pcitrif-method
setMethod("getRIFoutput", "pcitrif", function(x) x@step3$out)

#' @rdname pcitrif-class
#' @aliases getPCITinput_cond1,pcitrif-method
setMethod("getPCITinput_cond1", "pcitrif", function(x) x@step4$input_cond1)

#' @rdname pcitrif-class
#' @aliases getPCITinput_cond2,pcitrif-method
setMethod("getPCITinput_cond2", "pcitrif", function(x) x@step4$input_cond2)

#' @rdname pcitrif-class
#' @aliases getPCIToutput_cond1,pcitrif-method
setMethod("getPCIToutput_cond1", "pcitrif", function(x) x@step4$out_cond1)

#' @rdname pcitrif-class
#' @aliases getPCIToutput_cond2,pcitrif-method
setMethod("getPCIToutput_cond2", "pcitrif", function(x) x@step4$out_cond2)

#' @rdname pcitrif-class
#' @aliases getNetworkCond1,pcitrif-method
setMethod("getNetworkCond1", "pcitrif", function(x) x@step4$network_cond1)

#' @rdname pcitrif-class
#' @aliases getNetworkCond2,pcitrif-method
setMethod("getNetworkCond2", "pcitrif", function(x) x@step4$network_cond2)

#' @rdname pcitrif-class
#' @aliases getRIFkeyTF,pcitrif-method
setMethod("getRIFkeyTF", "pcitrif", function(x) x@step4$keytf)

#' @rdname pcitrif-class
#' @aliases getAnno,pcitrif-method
setMethod("getAnno", "pcitrif", function(x) x@step4$anno)
