// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// Centered distance matrix
arma::mat distance_center(const arma::vec& x) {
  int n = x.n_elem;
  arma::mat D(n, n);

  // Compute pairwise distances
  for (int i = 0; i < n; ++i)
    for (int j = i; j < n; ++j) {
      D(i, j) = std::abs(x(i) - x(j));
      D(j, i) = D(i,j);      
    }
      
  // Double centering
  arma::vec rowMeans = mean(D, 1);
  arma::vec colMeans = mean(D, 0).t();
  double grandMean = mean(rowMeans);

  for (int i = 0; i < n; ++i)
    for (int j = i; j < n; ++j) {
      D(i, j) = D(i, j) - rowMeans(i) - colMeans(j) + grandMean;
      D(j, i) = D(i, j);
    }
  return D;
}



//' Compute pairwise distance correlation metrics of each column to a vector
//' 
//' Fast computation of pairwise distance correlations.
//'
//' Note: To get the same result as from the energy package you need to take the square root of the results here.
//'
//' @param x A matrix. The number of rows should match the length of the vector y
//' @param y A vector
//' @return A vector with the same length as the number of columns in x. Each element contains the pairwise distance correlation to y.
//' @author Claus Ekstrom <claus@@rprimer.dk>
//' @examples
//' y <- rnorm(100)
//' x <- matrix(rnorm(100 * 10), ncol = 10)
//' pairwise_distance_correlation(x, y)
//'
//' @export
// [[Rcpp::export]]
arma::vec pairwise_distance_correlation(const arma::mat& x, const arma::vec& y) {
  int n = y.n_elem;
  int p = x.n_cols;
  arma::vec result(p);

  arma::mat A = distance_center(y);
  double dCovXX = accu(A % A) / (n * n);

  for (int j = 0; j < p; ++j) {
    arma::mat B = distance_center(x.col(j));
    double dCovYY = accu(B % B) / (n * n);
    double dCovXY = accu(A % B) / (n * n);

    result(j) = (dCovXX > 0 && dCovYY > 0) ? dCovXY / std::sqrt(dCovXX * dCovYY) : 0.0;
  }

  return result;
}
