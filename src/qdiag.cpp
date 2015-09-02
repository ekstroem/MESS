// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//' Fast extraction of matrix diagonal
//'
//' @description Fast extraction of matrix diagonal
//' @param x The matrix to extract the diagonal from
//' @return A vector with the diagonal elements
//' @details Note this function can only be used for extraction
//' @author Claus Ekstrom <claus@@rprimer.dk>
// [[Rcpp::export]]
NumericVector qdiag(NumericMatrix x) {

  int nrows = x.nrow();

  if (nrows != x.ncol())
    stop("The input matrix must be square");

  NumericVector res(nrows);

  for (int i=0; i<nrows; i++) {
    res(i) = x(i,i);
  }
 
  return res;
}

