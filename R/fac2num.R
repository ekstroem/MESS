#' Convert factor to numeric vector
#'
#' Converts the factor labels to a vector of numeric vectors
#'
#' Returns a vector of numeric values. Elements in the input factor that cannot be converted to numeric will produce NA.
#'
#' @param x A factor
#' @return Returns a numeric vector of the same length as x
#' @author Claus Ekstrom \email{claus@@rprimer.dk}
#' @keywords manip
#' @examples
#'
#' m1 <- matrix(rep(1, 16), 4)
#' m2 <- matrix(c(1, 0, .5, .5, 0, 1, .5, .5, .5, .5, 1, .5, .5, .5, .5, 1), 4)
#' m3 <- matrix(c(1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1), 4)
#' cmd(m1, m1)
#' cmd(m1, m2)
#' cmd(m2, m3)
#'
#' @export
fac2num <- function(x) {
    if (!is.factor(x))
        stop("input needs to be a factor")
    as.numeric(levels(x))[x]
}

