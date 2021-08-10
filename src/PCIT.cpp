// [[Rcpp::depends(RcppArmadillo)]]
//#include <Rcpp.h>
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace std;


//' @title Tolerance level between 3 pairwise correlations implemented in C/C++
//' 
//' @description Calculates the local tolerance for every trio of genes.
//' 
//' @param a Correlation value between the genes A and B.
//' @param b Correlation value between the genes B and C.
//' @param c Correlation value between the genes A and C.
//' @param tolType Calculation type for tolerance (1 for mean, 2 for min and 3 for max).
//'
//' @return Returns the value of tolerance.
//' 
//' @seealso See vignette for more details about the pairwise correlations.
//' 
//' @examples
//' tolerance(0.5, -0.65, 0.23, tolType = 1)
//' tolerance(0.5, -0.65, 0.23, tolType = 2)
//' tolerance(0.5, -0.65, 0.23, tolType = 3)
//' 
//' @rawNamespace useDynLib(CeTF)
//' @import Rcpp
//' @export
// [[Rcpp::export]]
float tolerance(float a, float b, float c, int tolType){
    float az, by, cx, tol;
    
    az = (a - b*c) / sqrt((1-pow(b, 2)) * (1-pow(c, 2)));
    by = (b - a*c) / sqrt((1-pow(a, 2)) * (1-pow(c, 2)));
    cx = (c - a*b) / sqrt((1-pow(a, 2)) * (1-pow(b, 2)));
    
    if (tolType == 1) {
        tol = (abs(az/a) + abs(by/b) + abs(cx/c)) / 3;
    }
    else if (tolType == 2) {
        tol = min(min(abs(az/a), abs(by/b)), abs(cx/c));
    }
    else if (tolType == 3) {
        tol = max(max(abs(az/a), abs(by/b)), abs(cx/c));
    } else {
        Rcpp::stop("invalid value of 'tolType'");
    }
    
    return tol;
}


//' @title A helper to calculate PCIT implemented in C/C++
//' 
//' @description Calculates the correlation matrix using PCIT algorithm
//'
//' @param cor A correlation matrix.
//' @param tolType Type of tolerance (default: 'mean') given the 3 pairwise correlations 
//' (see \code{\link{tolerance}}.
//' 
//' @return Correlation matrix resulted from PCIT algorithm.
//' 
//' @seealso (see \code{\link{PCIT}})
//' 
//' @examples
//' library(Matrix)
//' 
//' # loading a simulated normalized data
//' data('simNorm')
//' 
//' # calculating the correlation matrix
//' suppressWarnings(gene_corr <- cor(t(simNorm[1:30, ])))
//' gene_corr[is.na(gene_corr)] <- 0
//' 
//' # getting the PCIT correlation results for first 30 genes
//' results <- pcitC(cor = Matrix(gene_corr, sparse = TRUE), 
//'                 tolType = 1)
//' 
//' @rawNamespace useDynLib(CeTF)
//' @import Rcpp
//' @export
// [[Rcpp::export]]
arma::sp_mat pcitC(arma::sp_mat& cor, int tolType){
    float a, b, c, tol;
    float x, y, z;
    NumericVector xVals;
    
    arma::uword const nGenes = cor.n_rows;
    xVals = seq(0, nGenes-3);
    int nXVals = xVals.length();
    
    arma::sp_mat acor = cor;
    arma::sp_mat bcor = cor;
    
    for(int i=0; i<nXVals; i++) {
        x = xVals(i);
        for(y=x+1; y<nGenes-1; y++) {
            a = acor(x,y);
            for(z=y+1; z<nGenes; z++) {
                b = acor(x,z);
                c = acor(y,z);
                tol = tolerance(a, b, c, tolType);
                if(abs(a) < abs(b*tol) && abs(a) < abs(c*tol)) {
                    bcor(x,y) = 0;
                    bcor(y,x) = 0;
                }
                if(abs(b) < abs(a*tol) && abs(b) < abs(c*tol)) {
                    bcor(x,z) = 0;
                    bcor(z,x) = 0;
                }
                if(abs(c) < abs(a*tol) && abs(c) < abs(b*tol)) {
                    bcor(y,z) = 0;
                    bcor(z,y) = 0;
                }
            }
        }
    }
    return bcor;
}
