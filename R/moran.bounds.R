#' Function to compute extreme values of Moran's I 
#' 
#' This function computes the upper and lower bounds of Moran's I for a given spatial weighting matrix (stored in a \code{listw} object).
#' These values are obtained by the eigendecomposition of the spatial weighting matrix.
#'  
#' 
#' 
#' @param listw an object of class \code{listw}
#' @return A vector containing the maximum and minimum of Moran's I for a given spatial weighting matrix value returned
#' @author Stéphane Dray \email{stephane.dray@@univ-lyon1.fr}
#' @seealso \code{\link{mem}} \code{\link[spdep]{nb2listw}}
#' @references de Jong, P., Sprenger, C., & van Veen, F. (1984). On extreme values of Moran’s I and Geary's C. Geographical Analysis, 16(1), 17–24.
#' @keywords spatial
#' @examples
#' 
#' if(require("ade4", quietly = TRUE)){
#'  if(require("spdep", quietly = TRUE)){
#'      data(oribatid)
#'      nbtri <- tri2nb(as.matrix(oribatid$xy))
#'      lwB <- nb2listw(nbtri, style = "B")
#'      lwW <- nb2listw(nbtri, style = "W")
#'      scB <- mem(lwB)
#'      scW <- mem(lwW)

#'      moran.bounds(lwB)
#'      moran.mc(scB[,1], lwB, 9)
#'      moran.mc(scB[,69], lwB, 9)

#'      moran.bounds(lwW)
#'      moran.mc(scW[,1], lwW, 9)
#'      moran.mc(scW[,69], lwW, 9)
#'  }
#' }
#' 
#' @importFrom ade4 bicenter.wt
#' @importFrom spdep listw2mat
#' 
#' @export

moran.bounds <- function(listw){
    W <- listw2mat(listw)
    W <- (W + t(W)) / 2
    res <- nrow(W) / sum(W) *  range(eigen(bicenter.wt(W), only.values = TRUE))
    names(res) <- c("Imin", "Imax")
    return(res)
}

