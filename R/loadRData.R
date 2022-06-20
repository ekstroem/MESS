#' Load and extract object from RData file
#'
#' Loads and extracts an object from an RData file
#'
#' Returns an R object
#'
#' @param filename The path to the RData file
#' @return An R object
#' @author ricardo (from GitHub)
#' @seealso \code{\link{load}}
#' @keywords manip
#' @examples
#'
#' \dontrun{
#'   d <- loadRData("~/blah/ricardo.RData")
#' }
#'
#' @export loadRData
loadRData <- function(filename){
    #loads an RData file, and returns it
    load(filename)
    get(ls()[ls() != "filename"])
}
