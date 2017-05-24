#' Moment estimator
#'
#' To be filled out at a later date
#'
#' @param x input matrix, of dimension nobs x nvars; each row is an observation
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
#' @export
mestimate <- function(x, y) {
    
    res1 <- lm(y ~ x)             # Model 1
    res2 <- lm(I(y^2) ~ I(x^2))   # Model 2
    
    betahat <- (coef(res2)/coef(res1))[2]
    phat <- coef(res1)[2]/betahat

    list(phat=unname(phat),
         betahat=unname(betahat))

}


