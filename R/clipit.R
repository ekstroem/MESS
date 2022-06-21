#' Copy an object as R-code to the clipboard
#'
#' Copies an R object to the clipboard so it can be pasted in elsewhere.
#'
#' Returns nothing but will place the object in the clipboard
#'
#' @param x object to copy
#' @return Nothing but will put the R object into the clipboard as a side effect
#' @author Jonas LindeLøv posted on twitter \url{https://twitter.com/jonaslindeloev/status/1539182627554570240}. Copied shamelessly by Claus Ekstrom \email{claus@@rprimer.dk}
#' @keywords datagen
#' @examples
#'
#' clipit(mtcars$mpg)
#'
#' @importFrom clipr wite_clip
#' @export clipit
clipit <- function (x) 
{
    clipr::write_clip(capture.output(dput(x)))
}
