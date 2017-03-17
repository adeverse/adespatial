#' Construct a site-by-edge binary matrix
#' 
#' This function construct a site-by-edge binary matrix. It uses a set of sites
#' coordinates and a connexion diagram (object of class \code{nb} from the
#' \code{spdep} package). The 1s in the matrix represents the presence of a link
#' influencing a site, directly or indirectly, otherwise the values are 0s.
#' Graphically, the function is implemented such that the directional process is
#' considered to be going from the bottom to the top of the screen in the
#' graphical output of R. As such, the origin is underneath the set of points
#' representing the sites. Prepare the table of site coordinates accordingly.
#' 
#' The lowest site in the gradient is the one that will connect to the
#' fictitious site 0 to consider direction. Note that if there are multiple
#' lowest sites, they will all be connected to the fictitious site 0 to consider
#' direction.
#' 
#' @param nb.object Object of class \code{nb} from library spdep.
#' @param coords A three columns matrix or data frame. Colums 1: identifiers of
#'   the points (needs to be numeric). Columns 2 and 3: the X and Y coordinates
#'   of the points.
#' @param link A two columns matrix. Each row define an edge. Column 1: The site
#'   from which the edge starts. Column 2: the site to which the edge ends. All
#'   values in \code{link} need to be integers.
#' @param unit.angle Character. The measurement units in which the angle is
#'   defined: either "degrees" (default) or "radians".
#' @param rot.angle Numeric. Angle of the vector describing the process
#'   influencing the sites. This argument generate a rotation of the site
#'   coordinates. The set of coordinates is rotated counterclockwise. Negative
#'   values will produce a clockwise rotation.
#' @param rm.same.y Logical (\code{TRUE}, \code{FALSE}). Determines if the links
#'   perpendicular to the gradient should be removed. Default value:
#'   \code{TRUE}. If these links have already been removed this argument put to
#'   \code{TRUE} will make the function crash. See detail for more information.
#' @param plot.connexions Logical (\code{TRUE}, \code{FALSE}). Determines if the
#'   sites and the associated connexion diagram should be plotted after rotation
#'   of the coordinates by \code{gradient.angle}.
#'   
#' @details
#' 
#' The site-by-edge matrix created by this function and the list of edges
#' include the links associated to a fictitious site upstream of all other, see
#' Blanchet et al. (In press), for details. The decision regarding wether the
#' origin and the edges associated with it should be kept or removed is left to
#' the user. Removal of this site and of its associated edges can be done
#' manually after the construction of the site-by-edge matrix and of the list
#' edges. It can also be done when running the function \code{\link{aem}}.
#' 
#' If the connexion diagram was modified so that the links connecting sites that
#' are exactly perpendicular to the gradient have been removed or if there is no
#' sites exactly at the same level in the gradient, defining \code{rm.same.y} to
#' \code{TRUE} will generate an error.
#' 
#' If all the sites have the same y coordinates after rotation, e.g. a
#' horizontal transect perpendicular to the defined spatial asymmetry, this
#' analysis should not be used.
#' 
#' The argument \code{plot.connexions} will plot the sites (\code{coords}) in
#' black, after rotation, if any, and the connexion diagram (\code{nb.object}),
#' in red. The site labels are also plotted on the graph. To show the direction
#' of the spatial asymmetry considered by the function, a fictive site (in blue)
#' was added upstream. This fictive site is linked (blue edges) to the site(s)
#' that are the most upstream ones. Since this graph is generic, it might
#' sometimes look odd, however, the information given will remain the accurate.
#' 
#' @return
#' 
#' \item{\code{se.mat}}{A binary (n x k) matrix of site (n rows) by link edges
#' (k columns).} \item{\code{edges}}{A matrix describing the link edges. It has
#' 2 columns (from, to) and as many rows as there are edges. The edges linked to
#' the fictitious site of origin are found at the beginning of the list.}
#' 
#' @references
#' 
#' Blanchet F.G., P. Legendre and Borcard D. (2008) Modelling directional
#' spatial processes in ecological data. \emph{Ecological Modelling}, 215,
#' 325-336.
#' 
#' @author F. Guillaume Blanchet
#'   
#' @seealso \code{\link{aem}}
#'   
#' @importFrom spdep listw2sn
#' @importFrom graphics par
#' @importFrom graphics segments
#' @importFrom graphics text
#'   
#' @examples
#' 
#' ### Create an object of class nb (spdep)
#' if(require("spdep", quietly = TRUE)){
#' nb<-cell2nb(5,5,"queen")
#' 
#' ### Create fictitious geographical coordinates 
#' xy <- cbind(1:25,expand.grid(1:5,1:5))
#' 
#' ### Build a binary site-by-link matrix; remove the site which have identical Y coordinate
#' ### (by default argument: rm.same.y = TRUE)
#' bin.mat <- build.binary(nb,xy)
#' str(bin.mat)
#' 
#' ### Build a binary site-by-link matrix using the argument link: remove the site which
#' ### have identical Y coordinate (by default argument: rm.same.y = TRUE)
#' edges<-expand.grid(1,2:25)
#' bin.mat <- build.binary(coords=xy,link=edges)
#' str(bin.mat)
#' 
#' ### Build a binary site-by-link matrix, making the process affect the points at 
#' ### an angle of 45 degrees
#' bin.mat.45 <- build.binary(nb,xy, rot.angle=45)
#' str(bin.mat.45)
#' 
#' ### Build a binary site-by-link matrix, making the process affect the points at
#' ### an angle of pi/3 radians
#' bin.mat.pi3 <- build.binary(nb,xy,unit.angle="radians", rot.angle=pi/3)
#' str(bin.mat.pi3)
#' }
#' 
#' @keywords spatial
#' @export

`build.binary` <-
    function (nb.object=NULL, coords,link=NULL, unit.angle="degrees", rot.angle = 0, rm.same.y = TRUE, plot.connexions = TRUE) {
        if(is.character(coords[,1])){
            stop("The first column of 'coords' needs to be numeric")
        }else{
            coords<-as.matrix(coords)
        }
        if(!is.null(nb.object)){
            link <- listw2sn(nb2listw(nb.object, zero.policy = TRUE))[, 1:2]
        }else{
            if(is.null(link)){
                stop("both 'nb.object' and 'link' are NULL")
            }
        }
        link <- rm.double.link(link)
        n <- nrow(coords)
        if(missing(unit.angle)){
            unit.angle<-"degrees"
        }
        if (rot.angle == 0){
            if (unit.angle == "degrees" | unit.angle == "radians") {
                rot.angle<-0
            }else{
                stop("Units for angles be either 'degrees' or 'radians'")
            }
        }
        if (rot.angle != 0){
            if (unit.angle == "degrees") {
                rot.angle <- pi/180 * rot.angle
            }else{
                if (unit.angle == "radians") {
                    rot.angle <- rot.angle
                }else{
                    stop("Units for angles be either 'degrees' or 'radians'")
                }
            }
        }
        coords[, 2:3] <- round(rotation(as.matrix(coords[, 2:3]), rot.angle), digits = 8)
        if (rm.same.y == TRUE) {
            link <- remove.same.y(coords = coords, link = link)
        }
        
        low.y <- which(coords[, 3] == min(coords[, 3]))
        
        if (plot.connexions) {
            xy.range<-apply(coords[,2:3],2,range)
            
            xy.range.min<-xy.range[1,2]-((xy.range[2,2]-xy.range[1,2])/5)
            xy.range.max<-xy.range[2,2]
            
            par(mar=c(1,1,1,1))
            plot(coords[, 2:3], pch = 20, asp = 1, cex = 0.5, axes = FALSE, 
                xlab = "", ylab = "", ylim=c(xy.range.min,xy.range.max))
            segments(x0 = coords[link[, 1], 2], y0 = coords[link[, 
                1], 3], x1 = coords[link[, 2], 2], y1 = coords[link[, 
                    2], 3], col = "red")
            text(coords[,2:3],labels=as.character(coords[,1]),pos=2,col="red")
            
            site0<-c(mean(coords[,2]),xy.range.min)
            
            points(site0[1],site0[2],pch=19,col="blue")
            segments(x0=site0[1],y0=site0[2],x1=coords[low.y,2],y1=coords[low.y,3],col="blue")
            par(mar=c(5,4,4,2))
        }
        n.low.y <- length(low.y)
        link <- rbind(cbind(rep(0, n.low.y), low.y), as.matrix(link))
        nrow.link <- nrow(link)
        link <- cbind(1:nrow.link, link)
        points.order <- sort(coords[, 3], decreasing = FALSE, index.return = TRUE)$ix
        points.order <- c(0, points.order)
        link <- r.order.link(nrow.link, link, coords)
        link <- cbind(link, link[, 2])
        link.tmp <- link
        link2fac <- as.factor(link[, 2])
        matR <- mat.or.vec(1, n * nrow.link)
        mat <- .C("buildbinary", as.integer(nrow(link)), as.integer(link), 
            as.integer(points.order), as.integer(length(points.order)), 
            as.integer(n), matres = as.integer(matR), PACKAGE = "adespatial")$matres
        res.mat <- matrix(mat, nrow = n, ncol = nrow.link, byrow = TRUE)
        res.link <- link[, 2:3]
        colnames(res.link) <- c("from", "to")
        
        res<-list(res.mat, res.link)
        names(res)<-c("se.mat", "edges")
        return(res)
    }

