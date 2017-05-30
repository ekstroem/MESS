#' Moment estimator
#'
#' To be filled out at a later date
#'
#' @param x input vector, of dimension nobs x nvars; each row is an observation
#' vector.
#' @param y quantitative response variable of length nobs
#' @return Returns a list of 7 variables: \item{p.full }{The p-value for the
#' test of the full set of variables selected by the lasso (based on the OLS
#' estimates)} \item{ols.selected }{A vector of the indices of the non-zero
#' variables selected by \code{glmnet} sorted from (numerically) highest to
#' lowest based on their ols test statistic.} \item{p.maxols }{The p-value for
#' the maximum of the OLS test statistics} \item{lasso.selected }{A vector of
#' the indices of the non-zero variables selected by \code{glmnet} sorted from
#' (numerically) highest to lowest based on their absolute lasso coefficients.}
#' \item{p.maxlasso }{The p-value for the maximum of the lasso test statistics}
#' \item{lambda.orig }{The value of lambda used in the computations} \item{B
#' }{The number of permutations used}
#' @author Claus Ekstrom \email{ekstrom@@sund.ku.dk} 
#' @seealso \code{lm}
#' @keywords ~htests
#' @examples
#'
#' n <- 1000
#' p <- rbinom(n, size=1, prob=.20)
#' x <- rnorm(n)
#' y <- rnorm(n, mean=x)*p
#'
#' mestimate(x, y)
#' 
#' @export
mestimate <- function(x, y) {

    res1 <- coef(mfastLmCpp(y, cbind(x)))     # Model 1 (simple linear regression)
    res2 <- coef(mfastLmCpp(y^2, cbind(x^2))) # Model 2 (simple linear regression)

    ## Problem is coefficient for 1 is close to 0
    
    betahat <- res2/res1
    phat <- res1/betahat

    list(phat=unname(phat),
         betahat=unname(betahat))

}


