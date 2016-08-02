// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//' Fast marginal simple regresion analyses
//'
//' @description Fast computation of simple regression slopes for each predictor represented by a column in a matrix
//' @param y A vector of outcomes.
//' @param x A matrix of regressor variables. Must have the same number of rows as the length of y.
//' @param addintercept A logical that determines if the intercept should be included in all analyses (TRUE) or not (FALSE)
//' @return A data frame with two variables: coefficients and stderr that gives the slope estimate and corresponding standard error for each column in x.
//' @author Claus Ekstrom <claus@@rprimer.dk>
//' @examples
//' \dontrun{
//'   // Generate 100000 predictors and 100 observations
//'   x <- matrix(rnorm(100*100000))
//'   y <- rnorm(100, mean=x[,1])
//'   mfastLM_cpp(y, x)
//'
//' }
//' @export
// [[Rcpp::export]]
DataFrame mfastLm_cpp(NumericVector y, NumericMatrix x, bool addintercept=true) {
  arma::uword n = x.nrow(), k = x.ncol();
  int df = n-1;

  // Sanity checks
  if (y.size() != n) {
    stop("The length of y and the number of rows in x must match");
  }
  
  arma::mat X(x.begin(), n, k, false);
  arma::colvec Y(y.begin(), y.size(), false);
  arma::mat newX;
  arma::mat x0=arma::ones<arma::mat>(n,1);

  if (addintercept) {
    df = n-2;
  }

  arma::colvec rescoef = arma::zeros(k), resse = arma::zeros(k);
  arma::colvec coef, resid, stderrest;
  double sig2;
  
  for (arma::uword i=0; i<k; i++) {
    if (addintercept) {
      newX = join_rows(X.cols(i,i), x0);
    } else {
      newX = X.cols(i,i);
    }

    coef = arma::solve(newX, Y);
    rescoef(i) = coef(0);
    resid = Y - newX*coef;
    sig2 = arma::as_scalar(arma::trans(resid)*resid/df);
    stderrest = arma::sqrt(sig2 * arma::diagvec( arma::inv(arma::trans(newX)*newX)) );
    resse(i) = stderrest(0);
  }

  // create a new data frame and return it
  return DataFrame::create(Rcpp::Named("coefficients")=rescoef,
			   Rcpp::Named("stderr")=resse);
}

