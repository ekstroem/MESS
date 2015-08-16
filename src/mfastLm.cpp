// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//' Fast marginal simple regresion analyses
//'
//' @description Fast computation of the regression slopes for each predictor represented by a column in a matrix
//' @param y A vector of outcomes.
//' @param x A matrix of regressor variables. Must have the same number of rows as the length of y.
//' @param addintercept A logical that determines if the intercept should be included in all analyses (TRUE) or not (FALSE)
//' @return A data frame with two variables: coefficients and stderr that gives the slope estimate and corresponding standard error for each column in x.
//' @author Claus Ekstrøm <ekstrom@@sund.ku.dk>
//' @export
// [[Rcpp::export]]
DataFrame mfastLm_cpp(NumericVector y, NumericMatrix x, int addintercept) {
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


//' Fast marginal simple regresion analyses
//'
//' @description Fast computation of the regression slopes for each predictor represented by a column in a matrix
//' @return A data frame with two variables: coefficients and stderr that gives the slope estimate and corresponding standard error for each column in x.
//' @author Claus Ekstrøm <ekstrom@@sund.ku.dk>
//' @export
// [[Rcpp::export]]        
double lmm() {

  arma::uword N = 10000;
  arma::uword d = 5;

  arma::mat data(d, N, arma::fill::zeros);

  arma::vec mean0 = arma::linspace<arma::vec>(1,d,d);
  arma::vec mean1 = mean0 + 2;

  arma::uword i = 0;

  while(i < N)
    {
      if(i < N)  { data.col(i) = mean0 + arma::randn<arma::vec>(d); ++i; }
      if(i < N)  { data.col(i) = mean0 + arma::randn<arma::vec>(d); ++i; }
      if(i < N)  { data.col(i) = mean1 + arma::randn<arma::vec>(d); ++i; }
    }

  
  arma::gmm_diag model;
  model.learn(data, 2, arma::maha_dist, arma::random_subset, 10, 5, 1e-10, true);
  model.means.print("means:");

   double  scalar_likelihood = model.log_p( data.col(0)    );
  arma::rowvec     set_likelihood = model.log_p( data.cols(0,9));
  double overall_likelihood = model.avg_log_p(data);

  arma::uword   gaus_id  = model.assign( data.col(0),    arma::eucl_dist );
  arma::urowvec gaus_ids = model.assign( data.cols(0,9), arma::prob_dist );

  arma::urowvec hist1 = model.raw_hist (data, arma::prob_dist);
  arma::rowvec hist2 = model.norm_hist(data, arma::eucl_dist);

  /*
  model.save("my_model.gmm");
  */
  // the table is now initialized
  return 0;
}
