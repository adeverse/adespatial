#'@aliases aem.weight.time
#'  
#'@title Weight edges when constructing AEM variables
#'  
#'@description These functions construct a vector of weights that can be
#'  associated to the edges of the connexion diagram used as a basis to build
#'  AEM eigenfunctions. \code{aem.weight.edges} is general and can be used for 1 or
#'  2 dimensional problems. \code{aem.weight.time} is meant to be used only for time
#'  series. It is a wrapper for \code{aem.weight.edges}.
#'  
#'@param dates A vector of dates, class 'numeric' or 'Date'.
#'@param nb.object Object with class 'nb', computed by the spdep package,
#'  containing a list of neighbours for each sampling unit (site or time).
#'@param coords A three-column matrix or data frame. Column 1: identifiers of
#'  the points (must be numeric). Columns 2 and 3: the X and Y coordinates of
#'  the points.
#'@param distmat Class 'matrix' or 'dist' object containing a dissimilarity or
#'  distance matrix. (See details).
#'@param alpha Numeric. Exponent of the first weighting function. (See details).
#'@param beta Numeric. Exponent of the second weighting function. (See details).
#'@param max.d Numeric. Maximum distance for weighting. Default value if
#'  max.d=NULL: the maximum distance among a set of sites divided by 2 or the
#'  full span of a time series divided by 2 (not recommended in most problems,
#'  see details). A warning is given if \code{max.d = NULL} and the default
#'  value is used.
#'@param unit.angle Character. The measurement units in which the angle is
#'  defined: either "degrees" (default) or "radians".
#'@param rot.angle Numeric. Angle of the vector describing the process
#'  influencing the sites. This argument generates a rotation of the site
#'  coordinates. The set of coordinates is rotated counterclockwise. Negative
#'  values will produce a clockwise rotation.
#'@param rm.same.y Logical (\code{TRUE}, \code{FALSE}). Determines if the links
#'  perpendicular to the gradient should be removed. Default value: \code{TRUE}.
#'  If these links have already been removed, this argument put to \code{TRUE}
#'  will make the function crash. See details for more information.
#'@param plot.connexions Logical (\code{TRUE}, \code{FALSE}). Determines if the
#'  sites and the associated connexion diagram should be plotted after rotation
#'  of the coordinates by \code{gradient.angle}.
#'  
#'@details
#'
#'These functions should be used in close relationship with
#'\code{\link{aem.build.binary}}, consequently many of the arguments in this
#'function and in \code{\link{aem.build.binary}} are the same.
#'
#'The argument \code{distmat} may contain general forms of dissimilarity, for
#'example the difficulty of transfer of individuals, matter or energy among the
#'sampling units through space or time.
#'
#'In \code{aem.weight.edges}, two weighting functions, described in Legendre and
#'Legendre (2012, eqs. 114.3 and 14.4) have been implemented, where \eqn{d_{ij}}
#'is the distance between sites \eqn{i} and \eqn{j}:
#'
#'\tabular{ll}{ \code{Weighting function 1:} \tab \eqn{1 -
#'(d_{ij}/max(d))^\alpha}{1 - (d_{ij}/max(d))^alpha} \cr \code{Weighting
#'function 2:} \tab \eqn{1/d_{ij}^\beta}{1/d_{ij}^beta} }
#'
#'Also note that if a value is provided for \code{beta} (that is, if it is not
#'\code{NULL}), weighting function 2 is used regardless of whether \code{alpha}
#'is defined or not.
#'
#'In most applications, the default value of \code{max.d} is not optimal. A more
#'meaningful solution in many applications is to compute a Moran's I correlogram
#'(for univariate data) or a Mantel correlogram (for multivariate data), and
#'provide the distance where the correlation becomes 0 as the value for max.d.
#'
#'@return
#'
#'A vector of weights associating a value to each edge of the graph.
#'
#'@references
#'
#'Legendre, P. and L. Legendre (2012) \emph{Numerical Ecology}, 3rd English
#'edition. Elsevier Science BV, Amsterdam.
#'
#'Legendre, P. and O. Gauthier (2014) Statistical methods for temporal and
#'space-time analysis of community composition data. \emph{Proceedings of the
#'Royal Society B - Biological Sciences}, 281, 20132728.
#'
#'@author Olivier Gauthier, Pierre Legendre and F. Guillaume Blanchet
#'  
#'@seealso \code{\link{aem.build.binary}}, \code{\link[spdep]{sp.correlogram}},
#'  \code{\link[vegan]{mantel.correlog}}
#'  
#'  
#'@importFrom spdep listw2sn
#'@importFrom spdep cell2nb
#'@importFrom graphics par
#'@importFrom graphics segments
#'@importFrom graphics text
#'  
#' @examples
#'
#'### Time serie example
#'### Example - 12 dates (days from January 1st of year 1) 
#'### in a 6-year study starting September 5, 2000
#' if(require("spdep", quietly = TRUE)){
#' dates <- as.Date(c(129,269,500,631,864,976,1228,1352,1606,1730,1957,2087),origin="2000/1/1")
#'autocor.limit <- 522  # Limit of autcorrelation in the correlogram
#'
#'### Using aem.weight.time()
#'(wtime <- aem.weight.time(dates, alpha=2, max.d=autocor.limit))

#'### Using aem.weight.edges()
#'n <- length(dates)
#'nb <- cell2nb(1, n)
#'xy.dates <- cbind(1:n, rep(1, n), dates)
#'(wtime <- aem.weight.edges(nb, xy.dates, alpha=2, max.d=autocor.limit))
#'
#'n <- length(dates)
#'nb <- cell2nb(1, n)
#'xy.dates <- cbind(1:n, dates, rep(1, n)) ## Note the inversion of 'dates' and 'rep(1,n)'
#'wtime <- aem.weight.edges(nb, xy.dates, alpha=2, 
#'max.d=autocor.limit,rot.angle=90) # Note that 'rot.angle=90' was used
#'
#'### Spatial example using default d.max (notice the warning)
#'###########################################################################
#'nb<-cell2nb(5,5,"queen")
#'xy <- cbind(1:25,expand.grid(1:5,1:5))
#'(wspace <- aem.weight.edges(nb,xy))
#'}
#'
#' @keywords spatial
#' @keywords ts
#' @export

aem.weight.edges <-
    function(nb.object, coords, distmat=NULL, alpha=2, beta=NULL, max.d=NULL,
        unit.angle="degrees", rot.angle=0, rm.same.y=TRUE, plot.connexions=TRUE)
    {
        ### Olivier Gauthier, Pierre Legendre and F. Guillaume Blanchet - August 2013
        ##########################################################################################
        ### General check-up
        if(ncol(coords)!=3){
            stop("'coords' needs to have three columns")
        }
        if(is.character(coords[,1])){
            stop("The first column of 'coords' needs to be numeric")
        }else{
            coords<-as.matrix(coords)
        }
        
        if(rot.angle == 0){
            if (unit.angle == "degrees" | unit.angle == "radians") {
                rot.angle<-0
            }else{
                stop("Units for angles must be either 'degrees' or 'radians'")
            }
        }
        
        if(is.null(distmat)){
            links.length.mat <- as.matrix(dist(coords[,2:3]))
        }else{
            links.length.mat <- as.matrix(distmat)
        }
        if(is.null(max.d)){	
            max.d <- max(links.length.mat)/2
            print("The default value for 'max.d' was used. The weights may not be optimum")
        }
        
        ### Construct the links matrix
        links <- listw2sn(nb2listw(nb.object, zero.policy = TRUE))[, 1:2]
        links <- rm.double.link(links)
        
        ### Rotate the coordinates
        if(missing(unit.angle)){
            unit.angle<-"degrees"
        }
        
        if(rot.angle != 0){
            if (unit.angle == "degrees") {
                rot.angle <- pi/180 * rot.angle
            }else{
                if (unit.angle == "radians") {
                    rot.angle <- rot.angle
                }else{
                    stop("Units for angles must be either 'degrees' or 'radians'")
                }
            }
        }
        
        coords[, 2:3] <- round(rotation(as.matrix(coords[, 2:3]), rot.angle), digits = 8)
        
        ### Remove edges perpendicular to the directional gradient
        if(rm.same.y == TRUE) {
            links <- remove.same.y(coords = coords, link = links)
        }
        
        low.y <- which(coords[, 3] == min(coords[, 3]))
        
        ### Plot connexion diagram
        if(plot.connexions) {
            xy.range<-apply(coords[,2:3],2,range)
            
            xy.range.min<-xy.range[1,2]-((xy.range[2,2]-xy.range[1,2])/5)
            xy.range.max<-xy.range[2,2]
            
            par(mar=c(1,1,1,1))
            plot(coords[, 2:3], pch = 20, asp = 1, cex = 0.5, axes = FALSE, xlab = "", ylab = "", ylim=c(xy.range.min,xy.range.max))
            segments(x0 = coords[links[, 1], 2], y0 = coords[links[, 1], 3], x1 = coords[links[, 2], 2], y1 = coords[links[, 2], 3], col = "red")
            text(coords[,2:3],labels=as.character(coords[,1]),pos=2,col="red")
            
            site0<-c(mean(coords[,2]),xy.range.min)
            
            points(site0[1],site0[2],pch=19,col="blue")
            segments(x0=site0[1],y0=site0[2],x1=coords[low.y,2],y1=coords[low.y,3],col="blue")
            par(mar=c(5,4,4,2))
        }
        
        ### Reorganize the edges 
        n.low.y <- length(low.y)
        links <- rbind(cbind(rep(0, n.low.y), low.y), as.matrix(links))
        nrow.links <- nrow(links)
        links <- cbind(1:nrow.links, links)
        links <- r.order.link(nrow.links, links, coords)
        links <- links[, 2:3]
        
        ### Remove edges directly linked with site 0
        linkstorm<-unique(which(links==0,arr.ind=TRUE)[,1])
        links<-links[-linkstorm,]
        nrow.links <- nrow(links)
        
        ### Find the length of each edge
        links.length <- numeric(length=nrow.links)
        for(i in 1:nrow.links){
            links.length[i] <- links.length.mat[links[i,1],links[i,2]]
        }
        
        ### Calculate the weight associated to each edge
        if (is.null(beta)) {
            # First weighting function, Legendre & Legendre (2012, eq. 14.3)
            w <- 1 - (links.length/max.d)^alpha
            if(max(links.length) > max.d){
              stop(paste("'max.d' need to be larger than or equal to",max(links.length),"to have non-negative weights"))
            }
        } else {
            # Second weighting function, Legendre & Legendre (2012, eq. 14.4)
            w <- 1 / links.length^beta
        }
        w
    }
