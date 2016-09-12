#' Function to compute neighborhood based on the minimum spanning tree
#'
#' Compute mst and returns as a \code{nb} object
#'
#'
#' @param dxy A distance matrix based on spatial coordinates of samples
#' @return An object of class \code{nb}
#'
#' @author Stéphane Dray \email{stephane.dray@@univ-lyon1.fr}
#' @seealso  \code{\link[spdep]{graph2nb}}, \code{\link{give.thresh}}
#' @keywords spatial
#' @examples
#'
#' xy <- matrix(rnorm(60),30,2)
#' dxy <- dist(xy)
#' th <- give.thresh(dxy)
#' nb1 <- mst.nb(dxy)
#' nb1
#' wh1 <- which(as.matrix(dxy)==th,arr.ind=TRUE)
#' plot(nb1,xy,pch=20,cex=2,lty=3)
#' lines(xy[wh1[1,],1],xy[wh1[1,],2],lwd=2)
#' title(main="Maximum distance of the minimum spanning tree in bold")
#'
#' @importFrom ade4 mstree neig2nb
#' @export

"mst.nb" <-
    function(dxy) {
        mymst <- mstree(dxy)
        return(neig2nb(mymst))
    }
