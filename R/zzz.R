#' @importFrom utils packageDescription
.onAttach = function(libname, pkgname) {
    version = packageDescription(pkgname, fields = "Version")
    
    msg = paste0("========================================
", pkgname, 
        " version ", version, "
Bioconductor page: http://bioconductor.org/packages/CeTF/
Github page: https://github.com/cbiagii/CeTF or https://cbiagii.github.io/CeTF/
Documentation: http://bioconductor.org/packages/CeTF/
If you use it in published research, please cite:

Carlos Alberto Oliveira de Biagi Junior, Ricardo Perecin Nociti, Breno Osvaldo 
Funicheli, Patricia de Cassia Ruy, Joao Paulo Bianchi Ximenez, Wilson A Silva Jr.
CeTF: an R package to Coexpression for Transcription Factors using Regulatory 
Impact Factors (RIF) and Partial Correlation and Information (PCIT) analysis.
bioRxiv. 2020, DOI: 10.1101/2020.03.30.015784
========================================
")
    
    packageStartupMessage(msg)
}
