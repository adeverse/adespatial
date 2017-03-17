#' @describeIn weight.edges 
#' @export

weight.time <- function(dates, distmat=NULL,alpha=2, beta=NULL, max.d=NULL, unit.angle="degrees", rot.angle=0, rm.same.y=TRUE, plot.connexions=TRUE){
    #
    # dates : vector of dates.
    # alpha : exponent of the first weighting function
    # beta  : exponent of the second weighting function
    # max.d : maximum distance for weighting. Default value (not recommended in most 
    #         cases): the maximum distance among dates in the 'dates' vector. 
    #         A more meaningful solution in many applications is to compute a  
    #         Moran's I correlogram (for univariate data) or a Mantel correlogram 
    #         (for multivariate data), and use the distance where the correlation 
    #         becomes 0 as the value for max.d.
    #
    ### F. Guillaume Blanchet - August 2013
    ##########################################################################################
    n  <- length(dates)
    nb <- cell2nb(1,n)
    xy <- cbind(1:n, rep(1, n), dates)
    w  <- weight.edges(nb.object=nb, coords=xy, distmat=distmat, alpha=alpha, beta=beta, max.d=max.d, unit.angle=unit.angle, rot.angle=rot.angle, rm.same.y=rm.same.y, plot.connexions=plot.connexions)	
    w
}	
