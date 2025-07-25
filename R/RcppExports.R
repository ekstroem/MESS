# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Fast addition of vector to each row of matrix
#'
#' @description Fast addition of vector to each row of a matrix. This corresponds to t(t(x) + v)
#' @param x A matrix with dimensions n*k.
#' @param v A vector of length k.
#' @return A matrix of dimension n*k where v is added to each row of x
#' @author Claus Ekstrom <claus@@rprimer.dk>
#' @examples
#'
#' A <- matrix(1:12, ncol=3)
#' B <- c(1, 2, 3)
#'
#' add_torows(A, B)
#'
#' @export
add_torows <- function(x, v) {
    .Call('_MESS_add_torows', PACKAGE = 'MESS', x, v)
}

#' Fast binning of numeric vector into equidistant bins
#' 
#' Fast binning of numeric vector into equidistant bins
#'
#' Missing values (NA, Inf, NaN) are added at the end of the vector as the last bin returned if missinglast is set to TRUE
#' 
#' @param x A matrix of regressor variables. Must have the same number of rows as the length of y.
#' @param width The width of the bins
#' @param origin The starting point for the bins. Any number smaller than origin will be disregarded
#' @param missinglast Boolean. Should the missing observations be added as a separate element at the end of the returned count vector.
#' @return An list with elements counts (the frequencies), origin (the origin), width (the width), missing (the number of missings), and last_bin_is_missing (boolean) telling whether the missinglast is true or not.
#' @author Hadley Wickham (from SO: https://stackoverflow.com/questions/13661065/superimpose-histogram-fits-in-one-plot-ggplot) - adapted here by Claus Ekstrøm <claus@@rprimer.dk>
#' @examples
#'
#' set.seed(1)
#' x <- sample(10, 20, replace = TRUE)
#' bin(x, 15)
#' 
#' @export
bin <- function(x, width, origin = 0, missinglast = FALSE) {
    .Call('_MESS_bin', PACKAGE = 'MESS', x, width, origin, missinglast)
}

.chisq_test_cpp <- function(x, margin = 0L, statistic = 1L, B = 100000L) {
    .Call('_MESS_chisq_test_cpp', PACKAGE = 'MESS', x, margin, statistic, B)
}

#' Correlation matrix distance
#'
#' @description Computes the correlation matrix distance between two correlation matrices
#' @param x First correlation matrix
#' @param y Second correlation matrix
#' @return Returns the correlation matrix distance, which is a value between 0 and 1. The correlation matrix distance becomes
#' zero for equal correlation matrices and unity if they differ to a maximum extent.
#' @author Claus Ekstrom \email{claus@@rprimer.dk}
#' @references Herdin, M., and Czink, N., and Ozcelik, H., and Bonek, E. (2005). \emph{Correlation matrix distance, a meaningful measure for
#' evaluation of non-stationary mimo channels}. IEEE VTC.
#' @keywords univar
#' @examples
#'
#' m1 <- matrix(rep(1, 16), 4)
#' m2 <- matrix(c(1, 0, .5, .5, 0, 1, .5, .5, .5, .5, 1, .5, .5, .5, .5, 1), 4)
#' m3 <- matrix(c(1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1), 4)
#' cmd(m1, m1)
#' cmd(m1, m2)
#' cmd(m2, m3)
#'
#' @export cmd
cmd <- function(x, y) {
    .Call('_MESS_cmd', PACKAGE = 'MESS', x, y)
}

#' Apply cumsum to each column of matrix
#' 
#' Fast computation of apply(m, 2, cumsum)
#'
#' @param m A matrix
#' @return A matrix the same size as m with the column-wise cumulative sums.
#' @author Claus Ekstrom <claus@@rprimer.dk>
#' @examples
#'   # Generate a 100 by 10000 matrix
#'   x <- matrix(rnorm(100*10000), nrow=100)
#'   result <- colCumSum(x)
#'
#' @export
colCumSum <- function(m) {
    .Call('_MESS_colCumSum', PACKAGE = 'MESS', m)
}

#' Binning based on cumulative sum with reset above threshold
#' 
#' Fast binning of cumulative vector sum with new groups when the sum passes a threshold or the group size becomes too large
#'
#' Missing values (NA, Inf, NaN) are completely disregarded and pairwise complete cases are used f
#' 
#' @param x A matrix of regressor variables. Must have the same number of rows as the length of y.
#' @param threshold The value of the threshold that the cumulative group sum must not cross OR the threshold that each group sum must pass (when the argument cuwhatpassed is set to TRUE). 
#' @param cutwhenpassed A boolean. Should the threshold be the upper limit of the group sum (the default) or the value that each group sum needs to pass (when set to TRUE).
#' @param maxgroupsize An integer that defines the maximum number of elements in each group. NAs count as part of each group but do not add to the group sum. NULL (the default) corresponds to no group size limits.
#' @return An integer vector giving the group indices
#' @author Claus Ekstrom <claus@@rprimer.dk>
#' @examples
#'
#' set.seed(1)
#' x <- sample(10, 20, replace = TRUE)
#' cumsumbinning(x, 15)
#' cumsumbinning(x, 15, 3)
#' 
#' x <- c(3, 4, 5, 12, 1, 5, 3)
#' cumsumbinning(x, 10)
#' cumsumbinning(x, 10, cutwhenpassed=TRUE)
#'
#' @export
cumsumbinning <- function(x, threshold, cutwhenpassed = FALSE, maxgroupsize = NULL) {
    .Call('_MESS_cumsumbinning', PACKAGE = 'MESS', x, threshold, cutwhenpassed, maxgroupsize)
}

#' Fast distance covariance matrix
#'
#' @description Fast computation of the distance covariance between two matrices with the same number of rows.
#' @param x A matrix with dimensions n*k.
#' @param y A matrix with dimensions n*l.
#' @return A number representing the distance covariance between x and y
#' @author Claus Ekstrom <claus@@rprimer.dk>
#' @export
dCov <- function(x, y) {
    .Call('_MESS_dCov', PACKAGE = 'MESS', x, y)
}

#' Fast distance correlation matrix
#'
#' @description Fast computation of the distance correation matrix between two matrices with the same number of rows. Note that this is not the same as the correlation matrix distance that can be computed with the cmd function.
#' @param x A matrix with dimensions n*k.
#' @param y A matrix with dimensions n*l.
#' @return A number between 0 and 1 representing the distance covariance between x and y
#' @author Claus Ekstrom <claus@@rprimer.dk>
#' @export
dCor <- function(x, y) {
    .Call('_MESS_dCor', PACKAGE = 'MESS', x, y)
}

#' Compute pairwise distance correlation metrics of each column to a vector
#' 
#' Fast computation of pairwise distance correlations.
#'
#' Note: To get the same result as from the energy package you need to take the square root of the results here.
#'
#' @param x A matrix. The number of rows should match the length of the vector y
#' @param y A vector
#' @return A vector with the same length as the number of columns in x. Each element contains the pairwise distance correlation to y.
#' @author Claus Ekstrom <claus@@rprimer.dk>
#' @examples
#' y <- rnorm(100)
#' x <- matrix(rnorm(100 * 10), ncol = 10)
#' pairwise_distance_correlation(x, y)
#'
#' @export
pairwise_distance_correlation <- function(x, y) {
    .Call('_MESS_pairwise_distance_correlation', PACKAGE = 'MESS', x, y)
}

#' Fill down NA with the last observed observation
#'
#' @description Fill down missing values with the latest non-missing value
#' @param x A vector
#' @return A vector or list with the NA's replaced by the last observed value.
#' @author Claus Ekstrom <claus@@rprimer.dk>
#' @examples
#'
#' a <- c(1:5, "Howdy", NA, NA, 2:3, NA)
#' filldown(a)
#' filldown(c(NA, NA, NA, 3:5))
#'
#' @export
filldown <- function(x) {
    .Call('_MESS_filldown', PACKAGE = 'MESS', x)
}

#' Fast estimation of allele and genotype frequencies under Hardy-Weinberg equilibrium
#' 
#' Alleles are assumed to be numerated from 1 and up with no missing label. Thus if the largest value in either allele1 or allele2 is K then we assume that there can be at least K possible alleles.
#' Genotypes are sorted such the the smallest allele comes first, i.e., 2x1 -> 1x2, and 2x3 -> 2x3
#' 
#' @param allele1 An integer vector (starting with values 1 upwards) of first alleles
#' @param allele2 An integer vector (starting with values 1 upwards) of second alleles
#' @param min_alleles A minimum number of unique alleles available
#' @return A list with three variables: allele_freq for estimated allele frequencies, genotype_freq for estimated genotype_frequencies (under HWE assumption), obs_genotype is the frequency of the genotypes, available_genotypes is the number of available genotypes used for the estimation, and unique_alleles is the number of unique alleles (matches the length of allele_freq)
#' @author Claus Ekstrom <claus@@rprimer.dk>
#' @examples
#' al1 <- sample(1:5, size=1000, replace=TRUE, prob=c(.4, .2, .2, .1, .1))
#' al2 <- sample(1:5, size=1000, replace=TRUE, prob=c(.4, .2, .2, .1, .1))
#' hwe_frequencies(al1, al2)
#'
#' @export
hwe_frequencies <- function(allele1, allele2, min_alleles = 0L) {
    .Call('_MESS_hwe_frequencies', PACKAGE = 'MESS', allele1, allele2, min_alleles)
}

#' Kolmogorov-Smirnov goodness of fit test for cumulative discrete data
#'
#' The name of the function might change in the future so keep that in mind!
#'
#' @description Kolmogorov-Smirnov goodness of fit test for cumulative discrete data. 
#' @param x A vector representing the contingency table.
#' @param B The number of simulations used to compute the p-value.
#' @param prob A positive vector of the same length as x representing the distribution under the null hypothesis. It will be scaled to sum to 1. If NULL (the default) then a uniform distribution is assumed.
#' @details Simulation is done by random sampling from the null hypothesis.
#' @return A list of class "htest" giving the simulation results.
#' @author Claus Ekstrom <claus@@rprimer.dk>
#' @examples
#'
#' x <- 1:6
#' ks_cumtest(x)
#'
#' @export
ks_cumtest <- function(x, B = 10000L, prob = NULL) {
    .Call('_MESS_ks_cumtest', PACKAGE = 'MESS', x, B, prob)
}

#' Fast computation of maximum sum subarray
#'
#' @description Fast computation of the maximum subarray sum of a vector using Kadane's algorithm. The implementation handles purely negative numbers.
#' @param x A vector
#' @return A list with three elements: sum (the maximum subarray sum), start (the starting index of the subarray) and end (the ending index of the subarray)
#' @author Claus Ekstrom <claus@@rprimer.dk>
#' @examples
#'
#' maximum_subarray(1:4)
#' 
#' maximum_subarray(c(-2, 1, -3, 4, -1, 2, 1, -5, 4))
#'  
#' maximum_subarray(rnorm(100000))
#'
#' @export
maximum_subarray <- function(x) {
    .Call('_MESS_maximum_subarray', PACKAGE = 'MESS', x)
}

#' Fast marginal simple regresion analyses
#' 
#' Fast computation of simple regression slopes for each predictor represented by a column in a matrix
#'
#' No error checking is done
#' 
#' @param y A vector of outcomes.
#' @param x A matrix of regressor variables. Must have the same number of rows as the length of y. 
#' @param addintercept A logical that determines if the intercept should be included in all analyses (TRUE) or not (FALSE)
#' @return A data frame with three variables: coefficients, stderr, and tstat that gives the slope estimate, the corresponding standard error, and their ratio for each column in x.
#' @author Claus Ekstrom <claus@@rprimer.dk>
#' @examples
#' \dontrun{
#'   // Generate 100000 predictors and 100 observations
#'   x <- matrix(rnorm(100*100000), nrow=100)
#'   y <- rnorm(100, mean=x[,1])
#'   mfastLmCpp(y, x)
#'
#' }
#' @export
mfastLmCpp <- function(y, x, addintercept = TRUE) {
    .Call('_MESS_mfastLmCpp', PACKAGE = 'MESS', y, x, addintercept)
}

#' Compute Schur products (element-wise) of all pairwise combinations of columns in matrix
#'
#' Fast computation of all pairwise element-wise column products of a matrix.
#'
#' Note that the output order of columns corresponds to the order of the columns in x. First column 1 is multiplied with each of the other columns, then column 2 with the remaining columns etc. 
#'
#' @param x A matrix with dimensions r*c.
#' @param self A logical that determines whether a column should also be multiplied by itself.
#' @return A matrix with the same number of rows as x and a number of columns corresponding to c choose 2 (+ c if self is TRUE), where c is the number of columns of x. 
#' @author Claus Ekstrom <claus@@rprimer.dk>
#' @examples
#'
#' X <- cbind(rep(1, 4), 1:4, 4:1)
#' pairwise_Schur_product(X)
#' pairwise_Schur_product(X, self=TRUE)
#'
#' @export
pairwise_Schur_product <- function(x, self = FALSE) {
    .Call('_MESS_pairwise_Schur_product', PACKAGE = 'MESS', x, self)
}

#' Compute all pairwise combinations of indices
#'
#' Fast computation of indices of all pairwise element of a vector of length n.
#'
#' Note that the output order of columns corresponds to the order of the columns in x. First column 1 is multiplied with each of the other columns, then column 2 with the remaining columns etc. 
#'
#' @param n A number giving the number of elements to create all pairwise indices from
#' @param self A logical that determines whether a column should also be multiplied by itself.
#' @return A matrix with n*(n+1)/2 rows (if self=TRUE) or n*(n-1)/2 rows (if self=FALSE, the default) and two columns gicing all possible combinations of indices.
#' @author Claus Ekstrom <claus@@rprimer.dk>
#' @examples
#'
#' pairwise_combination_indices(3)
#' pairwise_combination_indices(4, self=TRUE)
#'
#' @export
pairwise_combination_indices <- function(n, self = FALSE) {
    .Call('_MESS_pairwise_combination_indices', PACKAGE = 'MESS', n, self)
}

#' Fast extraction of matrix diagonal
#'
#' @description Fast extraction of matrix diagonal
#' @param x The matrix to extract the diagonal from
#' @return A vector with the diagonal elements
#' @details Note this function can only be used for extraction
#' @author Claus Ekstrom <claus@@rprimer.dk>
#' @export qdiag
qdiag <- function(x) {
    .Call('_MESS_qdiag', PACKAGE = 'MESS', x)
}

#' Fast quadratic form computation
#'
#' @description Fast computation of a quadratic form  \eqn{t(x) * M * x}.
#' @param x A matrix with dimensions n*k.
#' @param M A matrix with dimenions n*n. If it is to be inverted then the matrix should be symmetric and positive difinite (no check is done for this)
#' @param invertM A logical. If set to TRUE then M will be inverted before computations (defaults to FALSE)
#' @param transposex A logical. Should the matrix be transposed before computations (defaults to FALSE).
#' @return A matrix with dimensions k * k giving the quadratic form
#' @author Claus Ekstrom <claus@@rprimer.dk>
#' @export
quadform <- function(x, M, invertM = FALSE, transposex = FALSE) {
    .Call('_MESS_quadform', PACKAGE = 'MESS', x, M, invertM, transposex)
}

#' Fast replication of a matrix
#'
#' @description Fast generation of a matrix by replicating a matrix row- and column-wise in a block-like fashion
#' @param x A matrix with dimensions r*c.
#' @param nrow An integer giving the number of times the matrix is replicated row-wise
#' @param ncol An integer giving the number of times the matrix is replicated column-wise
#' @return A matrix with dimensions (r*nrow) x (c*ncol)
#' @author Claus Ekstrom <claus@@rprimer.dk>
#' @examples
#'
#' m <- matrix(1:6, ncol=3)
#' repmat(m, 2)     # Stack two copies of m on top of each other
#' repmat(m, 2, 3)  # Replicate m with two copies on top and three copies side-by-side 
#'
#' @export
repmat <- function(x, nrow = 1L, ncol = 1L) {
    .Call('_MESS_repmat', PACKAGE = 'MESS', x, nrow, ncol)
}

#' Fast computation of trace of matrix product
#'
#' @description Fast computation of the trace of the matrix product trace(t(A) %*% B)
#' @param A A matrix with dimensions n*k.
#' @param B A matrix with dimenions n*k.
#' @return The trace of the matrix product
#' @author Claus Ekstrom <claus@@rprimer.dk>
#' @examples
#'
#' A <- matrix(1:12, ncol=3)
#' tracemp(A, A)
#'
#' @export
tracemp <- function(A, B) {
    .Call('_MESS_tracemp', PACKAGE = 'MESS', A, B)
}

