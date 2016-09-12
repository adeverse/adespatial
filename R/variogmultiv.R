#' Function to compute multivariate empirical variogram
#' 
#' Compute a multivariate empirical variogram. It is strictly equivalent to 
#' summing univariate variograms
#' 
#' 
#' @param Y A matrix with numeric data
#' @param xy A matrix with coordinates of samples
#' @param dmin The minimum distance value at which the variogram is computed 
#'   (i.e. lower bound of the first class)
#' @param dmax The maximum distance value at which the variogram is computed 
#'   (i.e. higher bound of the last class)
#' @param nclass Number of classes of distances
#' @return   A list: \item{d }{Distances (i.e. centers of distance classes).} 
#'   \item{var }{Empirical semi-variances.} \item{n.w }{Number of connections 
#'   between samples for a given distance.} \item{n.c }{Number of samples used 
#'   for the computation of the variogram.} \item{dclass }{Character vector with
#'   the names of the distance classes.}
#' @author Stéphane Dray \email{stephane.dray@@univ-lyon1.fr}
#' @references  Wagner H. H. (2003) Spatial covariance in plant communities: 
#'   integrating ordination, geostatistics, and variance testing. Ecology, 84, 
#'   1045--1057
#' @keywords spatial
#' @examples
#' 
#' if(require(ade4)){
#' data(oribatid)
#' # Hellinger transformation
#' fau <- sqrt(oribatid$fau / outer(apply(oribatid$fau, 1, sum), rep(1, ncol(oribatid$fau)), "*"))
#' # Removing linear effect
#' faudt <- resid(lm(as.matrix(fau) ~ as.matrix(oribatid$xy))) 
#' mvspec <- variogmultiv(faudt, oribatid$xy, nclass = 20)
#' mvspec
#' plot(mvspec$d, mvspec$var,type = 'b', pch = 20, xlab = "Distance", ylab = "C(distance)")
#' }
#' 
#' @export

"variogmultiv" <-
    function(Y,
        xy,
        dmin = 0,
        dmax = max(dist(xy)),
        nclass = 20) {
        dxy <- seq(dmin, dmax, le = nclass + 1)
        distgeo <- dist(xy)
        distfac <- cut(distgeo, bre = dxy)
        n.w <- table(distfac)
        colfac <-
            col(as.matrix(distgeo), as.factor = TRUE)[lower.tri(col(as.matrix(distgeo)))]
        rowfac <-
            row(as.matrix(distgeo), as.factor = TRUE)[lower.tri(row(as.matrix(distgeo)))]
        n.c <-
            apply(ifelse(table(distfac, colfac) + table(distfac, rowfac) > 0, 1, 0), 1, sum)
        
        Y2 <- dist(Y) ^ 2
        dxy2 <- dxy[-1]
        dxy <- dxy[-(nclass + 1)]
        
        res <- sapply(split(Y2, distfac), sum) / n.w / 2
        
        results <- list(
                d = (dxy + dxy2) / 2,
                var = as.vector(res),
                n.c = as.vector(n.c),
                n.w = as.vector(n.w),
                dclass = levels(distfac)
            )
        return(results)
    }
