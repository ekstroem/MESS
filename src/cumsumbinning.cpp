#include <Rcpp.h>
using namespace Rcpp ;

//' Binning based on cumulative sum with reset above threshold
//' 
//' Fast binning of cumulative vector sum with new groups when the sum passes a threshold or the group size becomes too large
//'
//' Missing values (NA, Inf, NaN) are completely disregarded and pairwise complete cases are used f
//' 
//' @param x A matrix of regressor variables. Must have the same number of rows as the length of y.
//' @param cutoff The value of the threshold that the cumulative group sum must not cross. 
//' @param maxgroupsize An integer that defines the maximum number of elements in each group. NULL (the default) corresponds to no group size.
//' @return An integer vector giving the group indices
//' @author Claus Ekstrom <claus@@rprimer.dk>
//' @examples
//'
//' set.seed(1)
//' x <- sample(10, 20, replace = TRUE)
//' cumsumbinning(x, 15)
//' cumsumbinning(x, 15, 3)
//' 
//' @export
// [[Rcpp::export]]
IntegerVector cumsumbinning(NumericVector x, double cutoff, Rcpp::Nullable<int> maxgroupsize=R_NilValue) {
  IntegerVector groupVec (x.size());
  int group = 1;
  int groupsize=0;
  double runSum = 0;


  int maxgroup;


  // Sanity checks
  if ((maxgroupsize.isNull())) {
    maxgroup=0;
  }
  else {     
    NumericVector xxx(maxgroupsize.get());
    maxgroup = xxx[0];
    if (maxgroup<1) {
      stop("maxgroupsize should be larger than 0");
    }
  }
  
  
  for (int i = 0; i < x.size(); i++) {
    if (!NumericVector::is_na(x[i]))
      runSum += x[i];
    groupsize++;

    if (runSum > cutoff) {
      group++;
      if (NumericVector::is_na(x[i])) {
	runSum = 0;
      } else {	
	runSum = x[i];
      }
      groupsize=0;
    }
    groupVec[i] = group;

    if ((maxgroup>0) && (groupsize==maxgroup)) {
      group++;
      runSum = 0;
      groupsize=0;      
    }
    
  }
  return groupVec;
}
