#' Compute the maximum distance of the minimum spanning tree based on a distance
#' matrix
#' 
#' It is used to select a truncation value for the dbMEM approach. It
#' returns the minimum value that keep all samples connected.
#' 
#' @param distxy A distance matrix (class \code{dist})
#'   
#' @return The maximum distance of the minimum spanning tree
#' @author Stéphane Dray \email{stephane.dray@@univ-lyon1.fr}
#' @export give.thresh
#' @importFrom ade4 mstree neig2mat
#'   
#' @examples
#' xy <- matrix(rnorm(60),30,2)
#' dxy <- dist(xy)
#' th <- give.thresh(dxy)

"give.thresh" <-
    function(distxy) {
        spanning <- mstree(distxy)
        return(max(neig2mat(spanning) * as.matrix(distxy)))
    }