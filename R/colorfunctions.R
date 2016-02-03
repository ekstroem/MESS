#' Shade an RGB color
#'
#' Shades an RBG color
#'
#' This function shades an RGB color and returns the shaded RGB color (with alpha channel added)
#'
#' @param x an RGB color
#' @param shade numeric value between 0 and 1. Zero means no change and 1 results in black
#' @author Claus Ekstrom \email{claus@@rprimer.dk}
#' @references Ekstrom, CT (2011) \emph{The R Primer}.
#' @keywords iplot
#' @examples
#'
#' newcol <- col.shade("blue")
#'
#' @export
col.shade <- function(x, shade=.5) {
    if (shade<0 | shade>1)
        stop("shade must be between 0 and 1")
    mat <- t(col2rgb(x, alpha=TRUE) * c(rep(1-shade, 3), 1))
    rgb(mat, alpha=mat[,4], maxColorValue=255)
}


#' Tint an RGB color
#'
#' Tints an RBG color
#'
#' This function tints an RGB color and returns the tinted RGB color (with alpha channel added)
#'
#' @param x an RGB color
#' @param tint numeric value between 0 and 1. Zero results in white and 0 means no change
#' @author Claus Ekstrom \email{claus@@rprimer.dk}
#' @references Ekstrom, CT (2011) \emph{The R Primer}.
#' @keywords iplot
#' @examples
#'
#' newcol <- col.tint("blue")
#'
#' @export
col.tint <- function(x, tint=.5) {
    if (tint<0 | tint>1)
        stop("shade must be between 0 and 1")
    mat <- t(col2rgb(x, alpha=TRUE)  +  c(rep(1-tint, 3), 0)*(255-col2rgb(x, alpha=TRUE)))
    rgb(mat, alpha=mat[,4], maxColorValue=255)
}
