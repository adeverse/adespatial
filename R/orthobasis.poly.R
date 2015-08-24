#' Function to compute polynomial of geographical coordinates
#' 
#' These functions compute polynomials of geographical coordinates. Polynomials 
#' functions are orthogonal and centred for the weights defined in \code{wt} 
#' (i.e., orthogonal to \code{wt}). It is the classical approach to perform 
#' trend surface analysis.
#' 
#' @param coords either a \code{Spatial*} object or a \code{matrix} with 
#'   geographic coordinates
#' @param degree the degree of the polynomial
#' @param wt a vector of weights. It is used to orthogonalize the polynomial 
#'   functions
#' @return an object of class \code{orthobasisSp} , subclass \code{orthobasis}
#' @author Stéphane Dray \email{stephane.dray@@univ-lyon1.fr}
#' @seealso \code{\link{mem}} \code{\link[ade4]{orthobasis}}
#' @references
#' 
#' Dray S., Pélissier R., Couteron P., Fortin M.J., Legendre P., Peres-Neto
#' P.R., Bellier E., Bivand R., Blanchet F.G., De Caceres M., Dufour A.B.,
#' Heegaard E., Jombart T., Munoz F., Oksanen J., Thioulouse J., Wagner H.H.
#' (2012). Community ecology in the age of multivariate multiscale spatial
#' analysis. \emph{Ecological Monographs} \bold{82}, 257--275.
#' 
#' @keywords spatial
#' @examples
#' 
#' if(require("ade4", quietly = TRUE)){
#' data(mafragh, package = "ade4")
#' pol2 <- orthobasis.poly(mafragh$Spatial) 
#' 
#' if(require("adegraphics", quietly = TRUE)){
#' plot(pol2, mafragh$Spatial)
#' }
#' }
#' 
#' 
#' @importMethodsFrom sp coordinates
#' @importFrom stats poly
#'   
#' @export
orthobasis.poly <- function(coords, degree = 2, wt = rep(1/nrow(coords), nrow(coords))){
    coords <- coordinates(coords) # coords could be a matrix, data.frame, SpatialPoints, SpatialPolygons, etc
    a0 <- poly(x = coords, degree = degree)
    poly.names <- colnames(a0)
    wt <- wt / sum(wt)
    a0 <- cbind(wt,a0)
    a0 <- qr.Q(qr(a0)) # to center the vectors
    a0 <- as.data.frame(a0[,-1]) / sqrt(wt)
    
    row.names(a0) <- row.names(coords)
    names(a0) <- sapply(strsplit(poly.names,"\\."), function(x) paste(c("X","Y"),x,sep ="", collapse="."))
    names(a0) <- gsub(".Y1","Y",gsub(".Y0","",gsub("X0.","",gsub("X1","X",names(a0)))))                
    attr(a0,"weights") <- wt
    attr(a0,"call") <- match.call()
    attr(a0,"class") <- c("orthobasisSp", "orthobasis","data.frame")
    
    return(a0)
}
