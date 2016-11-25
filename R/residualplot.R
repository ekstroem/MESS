#' @rdname residualplot
#' @export
residualplot.default <- function(x, y=NULL, candy=TRUE, bandwidth = 0.3, xlab="Fitted values", ylab="Std.res.", col.sd="blue", col.alpha=0.3,...) {

    if (is.null(y))
        stop("y must be specified")
  if (candy)
    plot(x, y, xlab = xlab, ylab = ylab,
         pch = 1 + 15 * (abs(y) > 1.96), ...)
  else plot(x, y, xlab = xlab, ylab=ylab, ...)
  if (candy) {

    # Set the colors
    if (col.alpha == 0)
      col.trans <- col.sd
    else col.trans <- sapply(col.sd, FUN = function(x) do.call(rgb,
                                     as.list(c(col2rgb(x)/255, col.alpha))))
    uniqx <- sort(unique(x))
    if (length(uniqx) > 3) {
      lines(smooth.spline(x, y, df = 3), lty = 2, lwd = 2,
            col = "black")
    }
    window <- bandwidth * (max(x) - min(x))/2
    vary <- length(uniqx)
    for (i in 1:length(uniqx)) {
      vary[i] <- 1.96 * sd(y[abs(x - uniqx[i]) <= window])
    }
    vary[is.na(vary)] <- 0
    polygon(c(uniqx, rev(uniqx)), c(vary, -(rev(vary))),
            col = col.trans, border = NA)
  }
  return(invisible(NULL))
}


#' @rdname residualplot
#' @export
residualplot.lm <- function(x, y, candy=TRUE, bandwidth = 0.3, xlab="Fitted values", ylab="Stud.res.", col.sd="blue", col.alpha=0.3,...) {
  y <- rstudent(x)
  x <- predict(x)
  residualplot(x, y, candy, bandwidth, xlab, ylab, col.sd, col.alpha, ...)
}





#' Plots a standardaized residual
#'
#' Plots a standardized residual plot from an lm object and provides additional
#' graphics to help evaluate the variance homogeneity and mean.
#'
#' Plots a standardized residual plot from an lm object and provides additional
#' graphics to help evaluate the variance homogeneity and mean.
#'
#' The blue area is a smoothed estimate of 1.96*SD of the standardized
#' residuals in a window around the predicted value. The blue area should
#' largely be rectangular if the standardized residuals have more or less the
#' same variance.
#'
#' The dashed line shows the smoothed mean of the standardized residuals and
#' should generally follow the horizontal line through (0,0).
#'
#' Solid circles correspond to standardized residuals outside the range from [-1.96; 1.96] while open circles are inside that interval. Roughly 5% of the observations should be outside the interval and the points should be evenly distributed.
#'
#' @aliases residualplot residualplot.lm residualplot.default
#' @param x lm object or a numeric vector
#' @param y numeric vector for the y axis values
#' @param candy logical. Should a lowess curve and local standard deviation of
#' the residual be added to the plot. Defaults to \code{TRUE}
#' @param bandwidth The width of the window used to calculate the local
#' smoothed version of the mean and the variance. Value should be between 0 and
#' 1 and determines the percentage of the window width used
#' @param xlab x axis label
#' @param ylab y axis label
#' @param col.sd color for the background residual deviation
#' @param col.alpha number between 0 and 1 determining the transprency of the
#' standard deviation plotting color
#' @param ... Other arguments passed to the plot function
#' @return Produces a standardized residual plot
#' @author Claus Ekstrom <claus@@rprimer.dk>
#' @seealso \code{\link{rstandard}}, \code{\link{predict}}
#' @keywords hplot
#' @examples
#'
#' # Linear regression example
#' data(trees)
#' model <- lm(Volume ~ Girth + Height, data=trees)
#' residualplot(model)
#'
#' @export
residualplot <- function(x, y=NULL, candy=TRUE, bandwidth = 0.3,
	                 xlab="Fitted values", ylab="Std.res.",
                         col.sd="blue", col.alpha=0.3,...) {
  UseMethod("residualplot")
}


