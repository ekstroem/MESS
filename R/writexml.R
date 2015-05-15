#' Write a data frame in XML format
#'
#' Writes the data frame to a file in the XML format.
#'
#' @param data the data frame object to save
#' @file the file name to be written to.
#' @na.rm logical. Should missing values be discarded from the written XML file (defaults to FALSE)
#' @return None
#' @details This function requires the \pkg{XML} package to be installed to function properly.
#'
#' @examples
#' data(trees)
#' write.xml(trees, file="mydata.xml")
#'
#' @author Claus Ekstrom, \email{claus@@rprimer.dk} based on previous work by Duncan Temple Lang.
#' @keyword file
#* @export write.xml

write.xml <- function(data, file=NULL, na.rm=FALSE) {
  if (!require(XML))
    stop("package XML must be installed")

  if(is.null(file))
    stop("filename not specified")

  if (!is.data.frame(data))
    stop("data must be a data frame")

  # Start empty XML document tree
  doc <- XML::newXMLDoc()          
  # Start by adding a document tag at the root of the XML file
  root <- XML::newXMLNode("document", doc=doc)
  
  # Make output invisible
  invisible(
    # Iterate over all rows
    lapply(1:nrow(data),                 
           function(rowi) {
             r <- XML::newXMLNode("row", parent=root)   # Create row tag
             for(var in names(data)) {   # Iterate over variables
               if (na.rm) {
	         if (! is.na(data[rowi, var])) {
	           XML::newXMLNode(var, data[rowi, var], parent = r)
                 }
               } else {
                 XML::newXMLNode(var, data[rowi, var], parent = r)
               }
             }
           }))            
  invisible(XML::saveXML(doc, file=file))
}
