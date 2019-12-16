#' The pcitrif Class
#'
#' The pcitrif class is data storage class that stores all the results of complete RIF and PCIT analysis.
#'
#' @slot step1 List of objects related to step 1 (Data Adjustment) of analysis.
#' @slot step2 List of objects related to step 2 (Differential Expression) of analysis.
#' @slot step3 List of objects related to step 3 (Regulatory Impact Factors) of analysis.
#' @slot step4 List of objects related to step 4 (Partial Correlation and Information Theory) of analysis.
#'
#' @importFrom methods new
#' @name pcitrif-class
#' @exportClass pcitrif
#'
pcitrif <- setClass(Class = "pcitrif", slots = list(step1 = "list", 
    step2 = "list", step3 = "list", step4 = "list"))



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

#' Method getGenes.
#' @name pcitrif-class
#' @rdname pcitrif-class
#' @exportMethod getGenes
setGeneric("getGenes", function(x) standardGeneric("getGenes"))

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

#' Method RIFinput.
#' @name pcitrif-class
#' @rdname pcitrif-class
#' @exportMethod RIFinput
setGeneric("RIFinput", function(x) standardGeneric("RIFinput"))

#' Method RIFout.
#' @name pcitrif-class
#' @rdname pcitrif-class
#' @exportMethod RIFout
setGeneric("RIFout", function(x) standardGeneric("RIFout"))

#' Method PCITinput1.
#' @name pcitrif-class
#' @rdname pcitrif-class
#' @exportMethod PCITinput1
setGeneric("PCITinput1", function(x) standardGeneric("PCITinput1"))

#' Method PCITinput2.
#' @name pcitrif-class
#' @rdname pcitrif-class
#' @exportMethod PCITinput2
setGeneric("PCITinput2", function(x) standardGeneric("PCITinput2"))

#' Method PCITout1.
#' @name pcitrif-class
#' @rdname pcitrif-class
#' @exportMethod PCITout1
setGeneric("PCITout1", function(x) standardGeneric("PCITout1"))

#' Method PCITout2.
#' @name pcitrif-class
#' @rdname pcitrif-class
#' @exportMethod PCITout2
setGeneric("PCITout2", function(x) standardGeneric("PCITout2"))

#' Method getNet1.
#' @name pcitrif-class
#' @rdname pcitrif-class
#' @exportMethod getNet1
setGeneric("getNet1", function(x) standardGeneric("getNet1"))

#' Method getNet2.
#' @name pcitrif-class
#' @rdname pcitrif-class
#' @exportMethod getNet2
setGeneric("getNet2", function(x) standardGeneric("getNet2"))

#' Method getKeyTF.
#' @name pcitrif-class
#' @rdname pcitrif-class
#' @exportMethod getKeyTF
setGeneric("getKeyTF", function(x) standardGeneric("getKeyTF"))

#' Method getAnno.
#' @name pcitrif-class
#' @rdname pcitrif-class
#' @exportMethod getAnno
setGeneric("getAnno", function(x) standardGeneric("getAnno"))



# setMethods
#' @rdname pcitrif-class
#' @return Return the raw data inputed
#' @aliases rawData,pcitrif-method
setMethod("rawData", "pcitrif", function(x) x@step1$raw)

#' @rdname pcitrif-class
#' @return Return the TPM data
#' @aliases tpmDat,pcitrif-method
setMethod("tpmData", "pcitrif", function(x) x@step1$tpm)

#' @rdname pcitrif-class
#' @return Return the filtered genes
#' @aliases getGenes,pcitrif-method
setMethod("getGenes", "pcitrif", function(x) x@step1$selected_genes)

#' @rdname pcitrif-class
#' @return Return the normalized data
#' @aliases normData,pcitrif-method
setMethod("normData", "pcitrif", function(x) x@step1$norm)

#' @rdname pcitrif-class
#' @return Return Differentially Expressed genes
#' @aliases getDE,pcitrif-method
setMethod("getDE", "pcitrif", function(x) x@step2$de_genes)

#' @rdname pcitrif-class
#' @return Return the TFs
#' @aliases getTF,pcitrif-method
setMethod("getTF", "pcitrif", function(x) x@step2$tf)

#' @rdname pcitrif-class
#' @return Return the RIF input
#' @aliases RIFinput,pcitrif-method
setMethod("RIFinput", "pcitrif", function(x) x@step3$input)

#' @rdname pcitrif-class
#' @return Return the RIF output
#' @aliases RIFout,pcitrif-method
setMethod("RIFout", "pcitrif", function(x) x@step3$out)

#' @rdname pcitrif-class
#' @return Return the PCIT input for condition1
#' @aliases PCITinput1,pcitrif-method
setMethod("PCITinput1", "pcitrif", function(x) x@step4$input_cond1)

#' @rdname pcitrif-class
#' @return Return the PCIT input for condition2
#' @aliases PCITinput2,pcitrif-method
setMethod("PCITinput2", "pcitrif", function(x) x@step4$input_cond2)

#' @rdname pcitrif-class
#' @return Return the PCIT output for condition1
#' @aliases PCITout1,pcitrif-method
setMethod("PCITout1", "pcitrif", function(x) x@step4$out_cond1)

#' @rdname pcitrif-class
#' @return Return the PCIT output for condition2
#' @aliases PCITout2,pcitrif-method
setMethod("PCITout2", "pcitrif", function(x) x@step4$out_cond2)

#' @rdname pcitrif-class
#' @return Return the network of interactions between gene-TFs for condition1
#' @aliases getNet1,pcitrif-method
setMethod("getNet1", "pcitrif", function(x) x@step4$network_cond1)

#' @rdname pcitrif-class
#' @return Return the network of interactions between gene-TFs for condition2
#' @aliases getNet2,pcitrif-method
setMethod("getNet2", "pcitrif", function(x) x@step4$network_cond2)

#' @rdname pcitrif-class
#' @return Return the main TFs after PCIT and RIF analysis
#' @aliases getKeyTF,pcitrif-method
setMethod("getKeyTF", "pcitrif", function(x) x@step4$keytf)

#' @rdname pcitrif-class
#' @return Return the annotation for each gene and TFs
#' @aliases getAnno,pcitrif-method
setMethod("getAnno", "pcitrif", function(x) x@step4$anno)
