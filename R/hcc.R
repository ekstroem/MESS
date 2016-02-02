#' Compute Hajnal-Curie-Cohen model
#'
#' XXXSplit a matrix into block diagonal sub matrices according to clusters and
#' XXXcombine the lower triangular parts into a vector
#'
#'
#' @param limit a square matrix
#' @param malefreq the frequency of male births (defaults to 0.51). Set to NULL to force equal male/female distribution among donor offspring
#' @param diag logical. Should the diagonal be included?
#' @return Returns a numeric vector containing the elements of the lower
#' triangular sub matrices.
#' @author Claus Ekstrom \email{claus@@ekstroem.dk}
#' @seealso \code{\link{lower.tri}}
#' @keywords manip
#' @examples
#'
#' m <- matrix(1:64, ncol=8)
#' cluster <- c(1, 1, 1, 1, 2, 2, 3, 3)
#' lower.tri.vector(m, cluster)
#'
#' @export lower.tri.vector

hcc <- function(limit = 12,
                Q = 6/7,
                malefreq = .51,
                familyfreq = 0,
                criteria = c("offspring", "families")) {

    ## Sanity check. The geography vectors should have the same length or be scalar

}
