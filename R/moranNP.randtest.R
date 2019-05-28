#' Function to compute positive and negative parts of Moran's index of spatial
#' autocorrelation
#' 
#' This function computes positive and negative parts of Moran's I statistic and
#' provide a testing procedure using random permutations. The functions compute
#' the Moran's eigenvector maps (MEM) and eigenvalues for the \code{listw}
#' object. If \code{alter = "greater"}, the statistic 'I+' is computed as the
#' sum of the products between positive eigenvalues and squared correlations
#' between \code{x} and associated MEMs. If \code{alter = "less"}, the statistic
#' 'I-' is computed as the sum of the products between negative eigenvalues and
#' squared correlations between \code{x} and associated MEMs. If \code{alter =
#' "two-sided"}, both statistics are computed.
#' 
#' @param x a \code{vector} with numeric data
#' @param listw an object of class \code{listw} created for example by
#'   \code{\link[spdep]{nb2listw}}
#' @param nrepet an integer indicating the number of permutations used in the 
#'   randomization procedure
#' @param alter a character string specifying the alternative hypothesis, must
#'   be one of "greater" (default), "less" or "two-sided"
#' @param \dots other arguments (e.g., \code{p.adjust.method}) to be passed to
#'   the code{\link[ade4]{as.krandtest}} function.
#' @return An object of class \code{randtest} (for unilateral test) or
#'   \code{krandtest} (for bilateral test)
#' @author Stéphane Dray \email{stephane.dray@@univ-lyon1.fr}
#' @seealso \code{\link{moran.randtest}}
#' @references Dray, S. (2011). A new perspective about Moran's coefficient:
#'   spatial autocorrelation as a linear regression problem. Geographical
#'   Analysis, 43, 127-141.
#' @keywords spatial
#' @examples
#' if(require("ade4", quietly = TRUE)  & require("spdep", quiet = TRUE)){
#' data(mafragh)
#' tests <- moranNP.randtest(mafragh$env[,1], nb2listw(mafragh$nb),
#'  alter = "two-sided", p.adjust.method = "holm")
#' tests
#' moran.randtest(mafragh$env[,1], nb2listw(mafragh$nb))$obs
#' sum(tests$obs)
#' }
#' 
#' @importFrom ade4 scalewt as.randtest
#' @importFrom spdep Szero
#' @useDynLib adespatial, .registration = TRUE 
#' @export moranNP.randtest

moranNP.randtest <- function(x, listw, nrepet = 999, alter = c("greater", "less", "two-sided"), ...) {
    
    alter <- match.arg(alter)
    z <- scalewt(x)
    scores <- mem(listw)
    res <- .C("testglobal", as.double(t(scores)), as.double(attr(scores, "values")), 
              as.integer(nrow(scores)), as.integer(ncol(scores)), as.double(z), 
              as.integer(nrepet), sim = as.double(rep(0,3*(nrepet+1))), PACKAGE = "adespatial")$sim
    
    res <- matrix(res, (nrepet + 1), 3, byrow = TRUE)
    colnames(res) <- c("I", "I+", "I-")
    res <- res / Szero(listw) / length(z)
    if(alter == "greater"){
        res <-  as.randtest(obs = res[1, 2], sim = res[-1, 2], alter = "greater", call = match.call())
    } else if(alter == "less"){
        res <-  as.randtest(obs = res[1, 3], sim = res[-1, 3], alter = "less", call = match.call())
    } else if(alter == "two-sided"){
        res <- as.krandtest(obs = res[1, 2:3], sim = res[-1, 2:3], alter = c("greater", "less"), call = match.call(), ...)
    }
    
    return(res)
}


