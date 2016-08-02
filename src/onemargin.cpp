// [[Rcpp::depends(RcppArmadillo)]]
// #include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;


// Compute chisq test statistic with no correction
double chisqtest(arma::mat x) {
  arma::colvec rowsum = sum(x, 1);
  arma::rowvec colsum = sum(x, 0);

  arma::mat expected = kron(colsum, rowsum)/sum(rowsum);

  return(  accu(square((x - expected))/expected));
}




//' Two-sided table test with fixed margins
//'
//' @description Test in a two-way contingency table with the row margin fixed. 
//' @param x A matrix representing the contingency table.
//' @param B The number of simulations used to compute the p-value.
//' @details Simulation is done by random sampling from the set of all tables with given row marginals, and works only if the marginals are strictly positive. Continuity correction is never used, and the statistic is quoted without it.
//' @return A list of class "htest" giving the simulation results.
//' @author Claus Ekstrom <claus@@rprimer.dk>
//' @examples
//'
//' m <- matrix(c(12, 4, 8, 6), 2)
//' chisq.test(m)
//' chisq.test(m, correct=FALSE)
//' fisher.test(m)
//' onemargintest(m)
//'
//' m2 <- matrix(c(9, 3, 3, 7), 2)
//' chisq.test(m2, simulate.p.value=TRUE)
//' fisher.test(m2)
//' onemargintest(m2)
//'
//' @export
// [[Rcpp::export]]
List onemargintest(NumericMatrix x, int B=10000) {

  arma::mat X(x.begin(), x.nrow(), x.ncol(), false);

  // Sanity checks
  if (B<1)
    Rcpp::stop("The number of permutation must be greater than 0");
  
  // Compute fixed row sizes
  arma::colvec n = sum(X, 1);

  // Compute column probabilities under H0
  arma::rowvec p = sum(X, 0)/sum(n);

  NumericVector PP = wrap(sum(X, 0)/sum(n));

  int nrows = n.size();
  int ncols = p.size();

  double originaltt = chisqtest(X);
  IntegerVector frame = seq_len(ncols) ;

  int larger=0;
  
  arma::mat testtable(nrows, ncols);
  
  for (int i=0; i<B; i++) {

    // Empty the simulated table
    testtable.zeros();    
    
    // Simulate each row    
    for (int j=0; j<nrows; j++) {
      IntegerVector res = RcppArmadillo::sample(frame, n(j), TRUE, PP) ;

      for (int k = 0; k < n(j); k++) {
	testtable(j, res(k)-1) +=1;
      }
    }

    /* Slightly faster but more convoluted

    int row=0;
    int idx = 0;
    for (int k = 0; k < N; k++) {
      if (idx == n(row)) {
	row +=1;
	idx=0;
      }
      testtable(row, res(k)-1) +=1;
      idx +=1;
    }


     */

    if (chisqtest(testtable) >= originaltt)
      larger +=1;

  }

  NumericVector statistic = NumericVector::create(_["X-squared"] = originaltt) ;

  Rcpp::List RVAL =  Rcpp::List::create(Rcpp::Named("method") = "Two-sided contingency table test with row margin fixed",
					Rcpp::Named("statistic") = statistic,
					Rcpp::Named("p.value") = (larger+1.0)/(B+1.0));

  RVAL.attr("class") = "htest";
  
  return(RVAL);
}

