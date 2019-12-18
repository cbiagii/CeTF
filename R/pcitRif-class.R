###################### pcitrif class

#' @title
#' The pcitrif Class
#'
#' @description
#' The pcitrif class is data storage class that stores all the
#' results from \code{\link{runAnalysis}} function.
#'
#' @slot step1 List of objects related to data adjustment.
#' @slot step2 List of objects related to differential expression analysis.
#' @slot step3 List of objects related to Regulatory Impact Factors analysis.
#' @slot step4 List of objects related to Partial Correlation and Information Theory analysis.
#'
#' @return Returns an pcitrif object.
#'
#' @importFrom methods new
#'
#' @name pcitrif-class
#'
#' @exportClass pcitrif
#'
pcitrif <- setClass(Class = "pcitrif", slots = list(step1 = "list",
    step2 = "list", step3 = "list", step4 = "list"))



###################### setGenerics
#' @title Accessors for the raw data in 'step1' slot of a
#' pcitrif object.
#' @description The step1 slot holds the data adjustment, that includes the
#' \code{raw data} (\code{\link{rawData}}), TPM data (\code{\link{tpmData}}),
#' filtered genes (\code{\link{getGenes}}) and normalized data
#' (\code{\link{normData}}).
#' @param x \code{\link{pcitrif-class}} object
#' @examples
#' # load the pcitrif class object resulted from runAnalysis function
#' data(pcitrifExample)
#'
#' rawData(pcitrifExample)
#' @seealso \code{\link{runAnalysis}}.
#' @name rawData
#' @rdname rawData-methods
#' @return Returns the raw data.
#' @exportMethod rawData
setGeneric("rawData", function(x) standardGeneric("rawData"))


#' @title Accessors for the TPM data in 'step1' slot of a
#' pcitrif object.
#' @description The step1 slot holds the data adjustment, that includes the
#' raw data (\code{\link{rawData}}), \code{TPM data} (\code{\link{tpmData}}),
#' filtered genes (\code{\link{getGenes}}) and normalized data
#' (\code{\link{normData}}).
#' @param x \code{\link{pcitrif-class}} object
#' @examples
#' # load the pcitrif class object resulted from runAnalysis function
#' data(pcitrifExample)
#'
#' tpmData(pcitrifExample)
#' @seealso \code{\link{countsToTPM}} and \code{\link{runAnalysis}}.
#' @name tpmData
#' @rdname tpmData-methods
#' @return Returns the TPM data.
#' @exportMethod tpmData
setGeneric("tpmData", function(x) standardGeneric("tpmData"))


#' @title Accessors for the filtered genes in 'step1' slot
#' of a pcitrif object.
#' @description The step1 slot holds the data adjustment, that includes the
#' raw data (\code{\link{rawData}}), TPM data (\code{\link{tpmData}}),
#' \code{filtered genes} (\code{\link{getGenes}}) and normalized data
#' (\code{\link{normData}}).
#' @param x \code{\link{pcitrif-class}} object
#' @examples
#' # load the pcitrif class object resulted from runAnalysis function
#' data(pcitrifExample)
#'
#' getGenes(pcitrifExample)
#' @seealso \code{\link{PCIT}} and \code{\link{runAnalysis}}.
#' @name getGenes
#' @rdname getGenes-methods
#' @return Returns the filtered genes.
#' @exportMethod getGenes
setGeneric("getGenes", function(x) standardGeneric("getGenes"))


#' @title Accessors for the normalized data in 'step1' slot of a pcitrif object.
#' @description The step1 slot holds the data adjustment, that includes the
#' raw data (\code{\link{rawData}}), TPM data (\code{\link{tpmData}}),
#' filtered genes (\code{\link{getGenes}}) and \code{normalized data}
#' (\code{\link{normData}}).
#' @param x \code{\link{pcitrif-class}} object
#' @examples
#' # load the pcitrif class object resulted from runAnalysis function
#' data(pcitrifExample)
#'
#' normData(pcitrifExample)
#' @seealso \code{\link{PCIT}} and \code{\link{runAnalysis}}.
#' @name normData
#' @rdname normData-methods
#' @return Returns the normalized data.
#' @exportMethod normData
setGeneric("normData", function(x) standardGeneric("normData"))


#' @title Accessors for the differentially expressed genes in 'step2' slot
#' of a pcitrif object.
#' @description The step2 slot holds the differential expression analysis,
#' that includes the \code{differentially expressed genes} (\code{\link{getDE}})
#' and the Transcription Factors (TFs) (\code{\link{getTF}}).
#' @param x \code{\link{pcitrif-class}} object
#' @examples
#' # load the pcitrif class object resulted from runAnalysis function
#' data(pcitrifExample)
#'
#' getDE(pcitrifExample)
#' @seealso \code{\link{expDiff}} and \code{\link{runAnalysis}}.
#' @name getDE
#' @rdname getDE-methods
#' @return Returns the Differentially Expressed (DE) genes.
#' @exportMethod getDE
setGeneric("getDE", function(x) standardGeneric("getDE"))


#' @title Accessors for the transcription factors in 'step2'
#' slot of a pcitrif object.
#' @description The step2 slot holds the differential expression analysis,
#' that includes the differentially expressed genes (\code{\link{getDE}})
#' and the \code{Transcription Factors (TFs)} (\code{\link{getTF}}).
#' @param x \code{\link{pcitrif-class}} object
#' @examples
#' # load the pcitrif class object resulted from runAnalysis function
#' data(pcitrifExample)
#'
#' getTF(pcitrifExample)
#' @seealso \code{\link{PCIT}} and \code{\link{runAnalysis}}.
#' @name getTF
#' @rdname getTF-methods
#' @return Returns the Transcription Factors (TFs).
#' @exportMethod getTF
setGeneric("getTF", function(x) standardGeneric("getTF"))


#' @title Accessors for the RIF analysis input in 'step3'
#' slot of a pcitrif object.
#' @description The step3 slot holds the Regulatory Impact Factors (RIF)
#' analysis, that includes the \code{input for RIF} (\code{\link{RIFinput}})
#' and output (\code{\link{RIFout}}).
#' @param x \code{\link{pcitrif-class}} object
#' @examples
#' # load the pcitrif class object resulted from runAnalysis function
#' data(pcitrifExample)
#'
#' RIFinput(pcitrifExample)
#' @seealso \code{\link{RIF}} and \code{\link{runAnalysis}}.
#' @name RIFinput
#' @rdname RIFinput-methods
#' @return Returns the input used for RIF analysis.
#' @exportMethod RIFinput
setGeneric("RIFinput", function(x) standardGeneric("RIFinput"))


#' @title Accessors for the RIF analysis output in 'step3'
#' slot of a pcitrif object.
#' @description The step3 slot holds the Regulatory Impact Factors (RIF)
#' analysis, that includes the input for RIF (\code{\link{RIFinput}})
#' and \code{output} (\code{\link{RIFout}}).
#' @param x \code{\link{pcitrif-class}} object
#' @examples
#' # load the pcitrif class object resulted from runAnalysis function
#' data(pcitrifExample)
#'
#' RIFout(pcitrifExample)
#' @seealso \code{\link{RIF}} and \code{\link{runAnalysis}}.
#' @name RIFout
#' @rdname RIFout-methods
#' @return Returns the output from RIF analysis.
#' @exportMethod RIFout
setGeneric("RIFout", function(x) standardGeneric("RIFout"))


#' @title Accessors for the PCIT input for condition1 in
#' 'step4' slot of a pcitrif object.
#' @description The step4 slot holds the Partial Correlation and Information
#' Theory (PCIT) analysis, that includes the PCIT input, output and the
#' network with gene-gene/gene-TFs interaction for condition1 and condition2,
#' besides the key TFs and the annotation (classification) of genes and TFs.
#' @param x \code{\link{pcitrif-class}} object
#' @examples
#' # load the pcitrif class object resulted from runAnalysis function
#' data(pcitrifExample)
#'
#' PCITinput1(pcitrifExample)
#' @seealso \code{\link{PCIT}} and \code{\link{runAnalysis}}.
#' @name PCITinput1
#' @rdname PCITinput1-methods
#' @return Returns the PCIT input for condition1.
#' @exportMethod PCITinput1
setGeneric("PCITinput1", function(x) standardGeneric("PCITinput1"))


#' @title Accessors for the PCIT input for condition2 in
#' 'step4' slot of a pcitrif object.
#' @description The step4 slot holds the Partial Correlation and Information
#' Theory (PCIT) analysis, that includes the PCIT input, output and the
#' network with gene-gene/gene-TFs interaction for condition1 and condition2,
#' besides the key TFs and the annotation (classification) of genes and TFs.
#' @param x \code{\link{pcitrif-class}} object
#' @examples
#' # load the pcitrif class object resulted from runAnalysis function
#' data(pcitrifExample)
#'
#' PCITinput2(pcitrifExample)
#' @seealso \code{\link{PCIT}} and \code{\link{runAnalysis}}.
#' @name PCITinput2
#' @rdname PCITinput2-methods
#' @return Returns the PCIT input for condition2.
#' @exportMethod PCITinput2
setGeneric("PCITinput2", function(x) standardGeneric("PCITinput2"))


#' @title Accessors for the PCIT output for condition1 in
#' 'step4' slot of a pcitrif object.
#' @description The step4 slot holds the Partial Correlation and Information
#' Theory (PCIT) analysis, that includes the PCIT input, output and the
#' network with gene-gene/gene-TFs interaction for condition1 and condition2,
#' besides the key TFs and the annotation (classification) of genes and TFs.
#' @param x \code{\link{pcitrif-class}} object
#' @examples
#' # load the pcitrif class object resulted from runAnalysis function
#' data(pcitrifExample)
#'
#' PCITout1(pcitrifExample)
#' @seealso \code{\link{PCIT}} and \code{\link{runAnalysis}}.
#' @name PCITout1
#' @rdname PCITout1-methods
#' @return Returns the PCIT output for condition1.
#' @exportMethod PCITout1
setGeneric("PCITout1", function(x) standardGeneric("PCITout1"))


#' @title Accessors for the PCIT output for condition2 in
#' 'step4' slot of a pcitrif object.
#' @description The step4 slot holds the Partial Correlation and Information
#' Theory (PCIT) analysis, that includes the PCIT input, output and the
#' network with gene-gene/gene-TFs interaction for condition1 and condition2,
#' besides the key TFs and the annotation (classification) of genes and TFs.
#' @param x \code{\link{pcitrif-class}} object
#' @examples
#' # load the pcitrif class object resulted from runAnalysis function
#' data(pcitrifExample)
#'
#' PCITout2(pcitrifExample)
#' @seealso \code{\link{PCIT}} and \code{\link{runAnalysis}}.
#' @name PCITout2
#' @rdname PCITout2-methods
#' @return Returns the PCIT output for condition2.
#' @exportMethod PCITout2
setGeneric("PCITout2", function(x) standardGeneric("PCITout2"))


#' @title Accessors for the network of interactions between key TFs and
#' genes for condition1 in 'step4' slot of a pcitrif object.
#' @description The step4 slot holds the Partial Correlation and Information
#' Theory (PCIT) analysis, that includes the PCIT input, output and the
#' network with gene-gene/gene-TFs interaction for condition1 and condition2,
#' besides the key TFs and the annotation (classification) of genes and TFs.
#' @param x \code{\link{pcitrif-class}} object
#' @examples
#' # load the pcitrif class object resulted from runAnalysis function
#' data(pcitrifExample)
#'
#' getNet1(pcitrifExample)
#' @seealso \code{\link{runAnalysis}}.
#' @name getNet1
#' @rdname getNet1-methods
#' @return Returns the network of interactions between gene and key TFs for condition1.
#' @exportMethod getNet1
setGeneric("getNet1", function(x) standardGeneric("getNet1"))


#' @title Accessors for the network of interactions between key TFs and
#' genes for condition2 in 'step4' slot of a pcitrif object.
#' @description The step4 slot holds the Partial Correlation and Information
#' Theory (PCIT) analysis, that includes the PCIT input, output and the
#' network with gene-gene/gene-TFs interaction for condition1 and condition2,
#' besides the key TFs and the annotation (classification) of genes and TFs.
#' @param x \code{\link{pcitrif-class}} object
#' @examples
#' # load the pcitrif class object resulted from runAnalysis function
#' data(pcitrifExample)
#'
#' getNet2(pcitrifExample)
#' @seealso \code{\link{runAnalysis}}.
#' @name getNet2
#' @rdname getNet2-methods
#' @return Returns the network of interactions between gene and key TFs for condition2.
#' @exportMethod getNet2
setGeneric("getNet2", function(x) standardGeneric("getNet2"))


#' @title Accessors for the key TFs in 'step4' slot of a pcitrif object.
#' @description The step4 slot holds the Partial Correlation and Information
#' Theory (PCIT) analysis, that includes the PCIT input, output and the
#' network with gene-gene/gene-TFs interaction for condition1 and condition2,
#' besides the key TFs and the annotation (classification) of genes and TFs.
#' @param x \code{\link{pcitrif-class}} object
#' @examples
#' # load the pcitrif class object resulted from runAnalysis function
#' data(pcitrifExample)
#'
#' getKeyTF(pcitrifExample)
#' @seealso \code{\link{runAnalysis}}.
#' @name getKeyTF
#' @rdname getKeyTF-methods
#' @return Returns the key TFs after PCIT and RIF analysis.
#' @exportMethod getKeyTF
setGeneric("getKeyTF", function(x) standardGeneric("getKeyTF"))


#' @title Accessors for the annotation of each gene and TFs in 'step4'
#' slot of a pcitrif object.
#' @description The step4 slot holds the Partial Correlation and Information
#' Theory (PCIT) analysis, that includes the PCIT input, output and the
#' network with gene-gene/gene-TFs interaction for condition1 and condition2,
#' besides the key TFs and the annotation (classification) of genes and TFs.
#' @param x \code{\link{pcitrif-class}} object
#' @examples
#' # load the pcitrif class object resulted from runAnalysis function
#' data(pcitrifExample)
#'
#' getAnno(pcitrifExample)
#' @seealso \code{\link{runAnalysis}}.
#' @name getAnno
#' @rdname getAnno-methods
#' @return Returns the annotation for each gene and TFs.
#' @exportMethod getAnno
setGeneric("getAnno", function(x) standardGeneric("getAnno"))





###################### setMethods
#' @rdname rawData-methods
#' @aliases rawData,pcitrif-method
setMethod("rawData", "pcitrif", function(x) x@step1$raw)

#' @rdname tpmData-methods
#' @aliases tpmData,pcitrif-method
setMethod("tpmData", "pcitrif", function(x) x@step1$tpm)

#' @rdname getGenes-methods
#' @aliases getGenes,pcitrif-method
setMethod("getGenes", "pcitrif", function(x) x@step1$selected_genes)

#' @rdname normData-methods
#' @aliases normData,pcitrif-method
setMethod("normData", "pcitrif", function(x) x@step1$norm)

#' @rdname getDE-methods
#' @aliases getDE,pcitrif-method
setMethod("getDE", "pcitrif", function(x) x@step2$de_genes)

#' @rdname getTF-methods
#' @aliases getTF,pcitrif-method
setMethod("getTF", "pcitrif", function(x) x@step2$tf)

#' @rdname RIFinput-methods
#' @aliases RIFinput,pcitrif-method
setMethod("RIFinput", "pcitrif", function(x) x@step3$input)

#' @rdname RIFout-methods
#' @aliases RIFout,pcitrif-method
setMethod("RIFout", "pcitrif", function(x) x@step3$out)

#' @rdname PCITinput1-methods
#' @aliases PCITinput1,pcitrif-method
setMethod("PCITinput1", "pcitrif", function(x) x@step4$input_cond1)

#' @rdname PCITinput2-methods
#' @aliases PCITinput2,pcitrif-method
setMethod("PCITinput2", "pcitrif", function(x) x@step4$input_cond2)

#' @rdname PCITout1-methods
#' @aliases PCITout1,pcitrif-method
setMethod("PCITout1", "pcitrif", function(x) x@step4$out_cond1)

#' @rdname PCITout2-methods
#' @aliases PCITout2,pcitrif-method
setMethod("PCITout2", "pcitrif", function(x) x@step4$out_cond2)

#' @rdname getNet1-methods
#' @aliases getNet1,pcitrif-method
setMethod("getNet1", "pcitrif", function(x) x@step4$network_cond1)

#' @rdname getNet2-methods
#' @aliases getNet2,pcitrif-method
setMethod("getNet2", "pcitrif", function(x) x@step4$network_cond2)

#' @rdname getKeyTF-methods
#' @aliases getKeyTF,pcitrif-method
setMethod("getKeyTF", "pcitrif", function(x) x@step4$keytf)

#' @rdname getAnno-methods
#' @aliases getAnno,pcitrif-method
setMethod("getAnno", "pcitrif", function(x) x@step4$anno)
