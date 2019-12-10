#' The pcitRif Class
#'
#' The pcitRif class is data storage class that stores all the results of complete RIF and PCIT analysis.
#'
#' @slot step1 List of objects related to step 1 (TPM Filter) of analysis.
#' @slot step2 List of objects related to step 2 (Differential Expression) of analysis.
#' @slot step3 List of objects related to step 3 (Regulatory Impact Factors analysis) of analysis.
#' @slot step4 List of objects related to step 4 (Partial Correlation and Information Theory analysis) of analysis.
#'
#' @importFrom methods new
#' @name pcitRif-class
#' @exportClass pcitRif
#'
pcitRif <- setClass(Class = "pcitRif", slots = list(step1 = "list", step2 = "list", step3 = "list", 
    step4 = "list"))
