#' Expand matrix or table to data frame
#'
#' Expands a contingency table to become a data frame where each observation becomes a combination of rows and columns.
#'
#' The procedure uses uniroot to find the root of a discontinuous function so
#' some errors may pop up due to the given setup that causes the root-finding
#' procedure to fail. Also, since exact binomial tests are used we have
#' discontinuities in the function that we use to find the root of but despite
#' this the function is usually quite stable.
#'
#' @param x A matrix
#' @return A data frame with the matrix expanded
#' @author Claus Ekstrom \email{claus@@rprimer.dk}
#' @keywords manip
#' @examples
#'
#' expand_matrix(diag(3))
#' m <- matrix(c(2, 1, 3, 0, 0, 2), 3)
#' expand_matrix(m)
#' @export
expand_matrix <- function(x) {

    if (! "matrix" %in% class(x))
        stop("needs matrix as input")

    x <- as.data.frame(as.table(x))
    x[rep(seq.int(1,nrow(x)), x$Freq), 1:2]
}
