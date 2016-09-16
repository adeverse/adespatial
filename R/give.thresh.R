#' Compute the maximum distance of the minimum spanning tree based on a distance
#' matrix
#' 
#' It is used to select a truncation value for the dbMEM approach. It
#' returns the minimum value that keep all samples connected.
#' 
#' @param matdist A distance matrix (class \code{dist} or \code{matrix})
#'   
#' @return The maximum distance in the minimum spanning tree.
#' @author Stéphane Dray \email{stephane.dray@@univ-lyon1.fr}
#' @export give.thresh
#' @importFrom ade4 mstree neig2mat
#'   
#' @examples
#' xy <- matrix(rnorm(60),30,2)
#' dxy <- dist(xy)
#' th <- give.thresh(dxy)

"give.thresh" <-
    function(matdist) {
        matdist <- as.dist(matdist)
        spanning <- mstree(matdist)
        return(max(neig2mat(spanning) * as.matrix(matdist)))
    }