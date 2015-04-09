#' Correlation matrix distance
#' 
#' Computes the correlation matrix distance between two correlation matrices
#' 
#' Returns a value between 0 and 1. The correlation matrix distance becomes
#' zero if the correlation matrices are equal up to a scaling factor and one if
#' they differ to a maximum extent.
#' 
#' @param m1 First correlation matrix
#' @param m2 Second correlation matrix
#' @return Returns the correlation matrix distance.
#' @author Claus Ekstrom \email{claus@@rprimer.dk}
#' @references Herdin, M., and Czink, N., and Ozcelik, H., and Bonek, E.
#' (2005). \emph{Correlation matrix distance, a meaningful measure for
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
cmd <- function(m1, m2) {
  if (!is.matrix(m1) || !is.matrix(m2))
    stop("requires two matrices to compute the cmd")

  if (any(dim(m1) - dim(m2) != 0))
    stop("the two matrices must have the same dimensions")
  
  1 - sum(diag(m1 %*% m2)) /(norm(m1, type="F")*norm(m2, type="F"))
}
