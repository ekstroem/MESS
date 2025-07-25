#' Copy an object as R-code to the clipboard
#'
#' Copies an R object to the clipboard so it can be pasted in elsewhere.
#'
#' Returns nothing but will place the object in the clipboard
#'
#' @param x object to copy
#' @return Nothing but will put the R object into the clipboard as a side effect
#' @author Jonas LindeLøv posted on twitter. Copied shamelessly by Claus Ekstrom \email{claus@@rprimer.dk}
#' @keywords datagen
#' @examples
#'
#'\dontrun{
#' clipit(mtcars$mpg)
#' }
#'
#' @importFrom clipr write_clip
#' @importFrom utils capture.output
#' @export clipit
clipit <- function (x) 
{
    clipr::write_clip(capture.output(dput(x)))
}
