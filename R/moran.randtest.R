#' Function to compute Moran's index of spatial autocorrelation
#' 
#' This function computes Moran's I statistic and provide a testing procedure using random permutations.
#' It is based on the \code{moran.mc} function of the \code{spdep} package. The \code{moran.randtest}
#' is slightly different as it allows to consider several variables (\code{x} can have more than one columns)
#' and its ouputs are objects of class \code{randtest} (one variable) or \code{krandtest} (several variables).
#'  
#' 
#' 
#' @param x a \code{vector}, \code{matrix} or \code{data.frame} with numeric data
#' @param listw an object of class \code{listw} created for example by \code{\link[spdep]{nb2listw}}
#' @param nrepet an integer indicating the number of permutations used in the 
#'  randomization procedure
#' @param \dots other arguments to be passed to the \code{\link[ade4]{as.randtest}} or code{\link[ade4]{as.krandtest}} functions. 
#' @return An object of class \code{randtest} (one variable) or \code{krandtest} (several variables)
#' 
#' @author Stéphane Dray \email{stephane.dray@@univ-lyon1.fr}
#' @seealso \code{\link[spdep]{moran.mc}}
#' @references Moran, P. A. P. (1950). Notes on continuous stochastic phenomena. Biometrika, 37, 17-23.
#' @keywords spatial
#' @examples
#' 
#' if(require("ade4", quietly = TRUE)  & require("spdep", quiet = TRUE)){
#' data(mafragh)
#' tests <- moran.randtest(mafragh$env, nb2listw(mafragh$nb))
#' tests
#' plot(tests)
#' 
#' }
#' 
#' @importFrom spdep moran.mc
#' @importFrom ade4 as.randtest as.krandtest
#' @export

moran.randtest <- function(x, listw, nrepet = 999, ... ){
    if(missing(listw) & inherits(x, "orthobasisSp"))
        if(!is.null(attr(x, "listw")))
            listw <- attr(x, "listw")
        
    x <- as.matrix(x)
    if(!is.numeric(x))
        stop("x should contain only numeric values")
    
    if(NCOL(x) == 1){
        res <- moran.mc(x, listw = listw, nsim = nrepet, zero.policy = TRUE)
        res <- as.randtest(obs = res$statistic, sim = res$res[-(nrepet + 1)], call = match.call(), ...)
    } else {
        res <- apply(x, 2, moran.mc, listw = listw, nsim = nrepet, zero.policy = TRUE)
        res <- as.krandtest(obs = sapply(res, function(x) x$statistic), sim = sapply(res, function(x) x$res[-(nrepet + 1)]), call = match.call(), ...)
    }
    
    return(res)
    
}