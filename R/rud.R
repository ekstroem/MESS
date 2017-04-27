#' Simulate randomized urn design
#'
#' Creates a plot of the crayon colors in \code{\link{brocolors}}
#'
#' @param n the number of individuals to randomize
#' @param alpha a non-negative integer vector of weights for each treatment group. The length of the vector corresponds to the number of treatment groups.
#' @param beta a non-negative integer    margin paramaters; vector of length 4 (see \code{\link[graphics]{par}})
#' @param labels margin paramaters; vector of length 4 (see \code{\link[graphics]{par}})
#' @param dataframe margin paramaters; vector of length 4 (see \code{\link[graphics]{par}})
#' @param startid margin paramaters; vector of length 4 (see \code{\link[graphics]{par}})
#'
#' @return None
#'
#' @examples
#' plot_crayons()
#'
#' @export
rud <- function(n, alpha=c(1, 1), beta=1, labels=seq(1, length(alpha)), dataframe=FALSE, startid=1) {

    ncat <- length(alpha)  # Number of categories
    
    if (ncat<2)
        stop("alpha mst have a length of at least 2")
    
    if (length(labels) != length(alpha))
        stop("length of labels must be the same as length of groups")
    
    urn <- alpha
    group <- seq(1, ncat)
    res <- integer(n)
    
    for (i in 1:n) {
        res[i] <- sample(group, 1, prob=urn/(sum(urn)))
        urn[-res[i]] <- urn[-res[i]] + beta
    }
    
    if (dataframe)
        data.frame(id=(1:n)+(startid-1), group=res, treatment=labels[res])
    else
        res
}
