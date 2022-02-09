#' convert cartesian coordinates to circular coordinates
#'
#' convert standard cartesian coordiates to coordinates in a circular plot
#'
#' @param x vector of x coordinates, or a 2-column matrix with columns x and y
#' @param y vector of y coordinates (ignored if x is a 2-column matrix)
#' @param xlim x-axis limits
#' @param ylim y-axis limits
#' @param rlim radius limits for plot
#'
#' @return 2-column matrix of (x,y) locations in circle plot
#'
#' @export

xy2circ <-
    function(x, y=NULL, xlim=NULL, ylim=NULL, rlim=c(1,2),
             start_angle=0)
{
    if(is.matrix(x)) {
        stopifnot(ncol(x) == 2)
        y <- x[,2]
        x <- x[,1]
    } else {
        stopifnot(length(x) == length(y))
    }

    if(is.null(xlim)) xlim <- range(x, na.rm=TRUE)
    if(is.null(ylim)) ylim <- range(y, na.rm=TRUE)

    # convert x's to angle in radians
    xnew <- (x - xlim[1])/(xlim[2]-xlim[1])*2*pi + start_angle

    # convert y's radius
    ynew <- (y - ylim[1])/(ylim[2]-ylim[1])*(rlim[2]-rlim[1]) + rlim[1]

    # convert angle, radius to x,y
    cbind(x=ynew*cos(xnew), y=ynew*sin(xnew))
}
