#' @importFrom utils packageDescription
.onAttach = function(libname, pkgname) {
    version = packageDescription(pkgname, fields = "Version")
    
    msg = paste0("========================================
", pkgname, 
        " version ", version, "
Bioconductor page: http://bioconductor.org/packages/CeTF/
Github page: https://github.com/cbiagii/CeTF
Documentation: http://bioconductor.org/packages/CeTF/

If you use it in published research, please cite:
Oliveira de Biagi Junior C, Perecin Nociti R, Osvaldo Funicheli B, 
Araujo da Silva Junior W (2020). CeTF: Coexpression for Transcription 
Factors using Regulatory Impact Factors and Partial Correlation and 
Information Theory analysis. R package version 0.99.13.
========================================
")
    
    packageStartupMessage(msg)
}
