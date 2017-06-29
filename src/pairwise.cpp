// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//' Compute Schur products (element-wise) of all pairwise combinations of columns in matrix
//'
//' Fast computation of all pairwise element-wise column products of a matrix.
//'
//' Note that the output order of columns corresponds to the order of the columns in x. First column 1 is multiplied with each of the other columns, then column 2 with the remaining columns etc. 
//'
//' @param x A matrix with dimensions r*c.
//' @param self A logical that determines whether a column should also be multiplied by itself.
//' @return A matrix with the same number of rows as x and a number of columns corresponding to c choose 2 (+ c if self is TRUE), where c is the number of columns of x. 
//' @author Claus Ekstrom <claus@@rprimer.dk>
//' @examples
//'
//' X <- cbind(rep(1, 4), 1:4, 4:1)
//' pairwise_Schur_product(X)
//' pairwise_Schur_product(X, self=TRUE)
//'
//' @export
// [[Rcpp::export]]
NumericMatrix pairwise_Schur_product(NumericMatrix x, bool self=false) {

  arma::mat X(x.begin(), x.nrow(), x.ncol(), false);
  arma::mat res(x.nrow(), (x.ncol()*(x.ncol()-1)/2  + ((self) ? x.ncol() : 0 )));

  int index = 0;
  for (int i=0; i<(x.ncol()-1 + ((self) ? 1 : 0)); i++) {
    for (int j=(i+1 + ((self) ? -1 : 0)); j<x.ncol(); j++) {
      res.col(index) = X.col(i) % X.col(j);
      index++;
    }
  }
  
  return wrap(res);
}

