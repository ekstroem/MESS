#' Screen variable before penalized regression
#'
#' Expands a contingency table to a data frame where each observation in the table becomes a single observation in the data frame with corresponding information for each for each combination of the table dimensions.
#'
#' Note that no standardization is done (not necessary?)
#'
#' @param x A table or matrix
#' @param y A vector out outcomes
#' @param lambda a vector of values penalized ....
#' @param method a string giving the method used for screening. Currently not used
#' @return A list with two elements: lambda which
#' @author Claus Ekstrom \email{claus@@rprimer.dk}
#' @keywords manip
#' @examples
#'
#' x <- matrix(rnorm(50*100), nrow=50)
#' y <- rnorm(50)
#' screen_variables(x, y, lambda=c(0, 1, 2))
#'
#' @export
screen_variables <- function(x, y, lambda=0, method=c("xxx")) {

    ## Sanity checks
    if (!any(c("matrix", "table") %in% class(x)))
        stop("needs matrix or table as input")

    ##    ndim <- length(dim(x))

    res <- abs(as.vector(crossprod(x, y)))
    lambdamax <- max(res)
    # discard <- res<(2*lambda-lambdamax)

    discard <- outer(res, 2*lambda-lambdamax, "<")

    # The code below might be modified and fast but requires data.table

    list(lambda=lambda, selected=unname(apply(discard, 2, which)))
}
