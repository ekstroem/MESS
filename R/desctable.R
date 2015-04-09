desctable.factor <- function(object, group=NULL, proportions=FALSE, useNA="ifany", ...) {
    if ((length(group) == length(object)) & (!is.null(group))) {
        group <- factor(group)
        u <- table(object, factor(group))
        p <- chisq.test(u)$p.value
    } else {
        u <- table(object, useNA = useNA)
        p <- NULL
    }
    if (length(u) == 0 || sum(u) == 0) { res <- "NA" }
    else { 
      if (proportions) {
        res <- prop.table(u)
      } else {
        res <- u
      }
    }
    list(output=res, p=p)
}

desctable.character <- function(object, group=NULL, ...) {
    desctable.factor(factor(object), group=group)
}


desctable.numeric <- function(object, group=NULL, proportions=FALSE, ...) {
    if (is.null(group)) {        
        mu <- mean(object, na.rm=TRUE)
        su <- sd(object, na.rm=TRUE)
        if (is.na(mu)) { res <- "NA" }
        else { res <- sprintf( ifelse( abs(mu) < 1, "%0.2f (%0.2f; -%0.2f)", "%0.2f (%0.2f; %0.2f)" ), mu, mu - 1.96*su, mu + 1.96*su ) }
        p <- NULL
    } else {
        group <- factor(group)
        mus <- by(object, group, mean, na.rm=TRUE)
        sds <- by(object, group, sd, na.rm=TRUE)
        res <- sapply(1:length(mus), function(xxx) { sprintf( "%0.2f (%0.2f; %0.2f)", mus[xxx], mus[xxx] - 1.96*sds[xxx], mus[xxx] + 1.96*sds[xxx] )} )
     #   names(res) <- unique(group)
        p <- oneway.test(object ~ group)$p.value
    }

#    res
    list(output=res, p=p)
}

desctable.integer <- function(object, group=NULL, proportions=FALSE, ...) {
    desctable.numeric(object)
}

# Same function, for data.frames
desctable.data.frame <- function(object, group=NULL, proportions=FALSE, ...) {

#    if (is.null(group)) {    
        res <- lapply(object, function(x) { desctable(x, group=group)})
#    } else {
#        res <- t(aggregate(object, list(group), desctable, simplify=TRUE))
#        res <- t(aggregate(object, list(group), desctable, simplify=TRUE))
        
#    }
#    data.frame(res)
        res
}

desctable.default <- function(object, group=NULL, proportions=FALSE, ...) {
    warning(paste("Cannot make description for class", class(object)))
    
}




#' Create a descriptive table of summaries
#' 
#' Creates a descriptive table of summary statistics
#' 
#' 
#' @aliases desctable desctable.character desctable.factor desctable.default
#' desctable.data.frame desctable.integer desctable.numeric
#' @param object an object for which a descriptive summary table is desired
#' @param group A factor the same length as object (or the elements of object
#' if that is a data frame) that is used to split up the summary for each level
#' of group. If set to NULL (the default) then every observation is assumed to
#' be in one group
#' @param proportions Logical. Should proportions be listed for categorical
#' data (defaults to FALSE)
#' @param \dots any other arguments passed of functions
#' @return A list with length the same as object with items \item{output
#' }{Output for use with the descriptive table} \item{p }{p-value for testing
#' equality across groups (set to NULL if groups=NULL)}
#' @author Claus Ekstrom \email{claus@@rprimer.dk}
#' @seealso \code{\link{summary}}
#' @keywords print
#' @examples
#' 
#' desctable(ToothGrowth, ToothGrowth$supp)
#' 
#' @export desctable
desctable <- function(object, group=NULL, proportions=FALSE, ...) {

    UseMethod("desctable")
   
}

