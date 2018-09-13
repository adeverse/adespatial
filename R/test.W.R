#' Function to compute and test eigenvectors of spatial weighting matrices
#' 
#' This function is now deprecated. Please try the new \code{\link{listw.candidates}} and 
#' \code{\link{listw.select}} functions.
#' 
#' This function is a user-friendly way to compute and test eigenvectors for
#' various definitions of spatial weighting matrices. It combines calls to the 
#' functions \code{scores.listw} and \code{ortho.AIC}. It allows to test various
#' definitions of the spatial weighting matrix and return results of 
#' \code{scores.listw} for the best one.
#' 
#' @details This functions allows to test one binary spatial weighting matrix 
#'   (if only Y and nb are provided). It allows also to test a weighting 
#'   function based on distances (if f is provided) and a weighting function 
#'   with different values of parameters if other arguments of \code{f} are 
#'   provided.
#'   
#' @param Y A matrix with response variables (univariate or multivariate 
#'   response)
#' @param nb An object of the class \code{nb} created by functions of the 
#'   \code{spdep} package
#' @param xy Coordinates of the samples, this argument is optional and is 
#'   required only if the argument \code{f} is not null.
#' @param MEM.autocor A string indicating if all MEM must be returned or only 
#'   those corresponding to positive or negative autocorrelation
#' @param f A function of the distance that can be used as a weighting spatial 
#'   function. This argument is optional
#' @param \dots Others arguments for the function \code{f}. It defines the range
#'   of parameters which will be tested
#' @return A list with the following elements: \item{all }{A data.frame where 
#'   each row correspond to one spatial weighting matrix tested. It contains 
#'   value of parameteres tested and corrected AIC and number of orthogonal 
#'   vectors for the best model.} \item{best }{A list containing results of 
#'   scores.listw and ortho.AIC of the best spatial weighting matrix according 
#'   to corrected AIC.}
#' @author Stéphane Dray \email{stephane.dray@@univ-lyon1.fr}
#' @seealso   \code{\link{ortho.AIC}}, \code{\link{scores.listw}}
#' @references Dray S., Legendre P. and Peres-Neto P. R. (2006) Spatial 
#'   modeling: a comprehensive framework for principal coordinate analysis of 
#'   neighbor matrices (PCNM). Ecological Modelling, 196, 483--493
#' @keywords spatial
#' @examples
#' 
#' if(require(ade4) & require(spdep)){
#' 
#' data(oribatid)
#' # Hellinger transformation
#' fau <- sqrt(oribatid$fau / outer(apply(oribatid$fau, 1, sum), rep(1, ncol(oribatid$fau)), "*"))
#' # remove gradient effect
#' faudt <- resid(lm(as.matrix(fau) ~ as.matrix(oribatid$xy)))
#' 
#' # test a binary spatial weighting matrix
#' nbtri <- tri2nb(as.matrix(oribatid$xy))
#' tri.res <- test.W(faudt, nbtri)
#' 
#' maxi <- max(unlist(nbdists(nbtri, as.matrix(oribatid$xy))))
#' 
#' # test a simple spatial weighting function of the distance
#' f1 <- function(x) {1-(x)/(maxi)}
#' tri.f1 <- test.W(faudt, nbtri, f = f1, xy = as.matrix(oribatid$xy))
#' 
#' # test a spatial weighting function with various values of parameters
#' f2 <- function(x,dmax,y) {1-(x^y)/(dmax)^y}
#' tri.f2 <- test.W(faudt,nbtri, f = f2, y = 2:10, dmax = maxi, xy = as.matrix(oribatid$xy))
#' }
#' 
#' @importFrom spdep nb2listw
#' @export
#' 
"test.W" <-
    function(Y,
        nb,
        xy,
        MEM.autocor = c("all", "positive", "negative"),
        f = NULL,
        ...) {
        .Deprecated(new = "listw.select", package = "adespatial", 
            msg = "This function is now deprecated. Please try the new 'listw.select' function.")
        
        mycall <- pairlist(...)
        res <- list()
        MEM.autocor <- match.arg(MEM.autocor)
        if (!(is.null(f))) {
            nbdist <- nbdists(nb, as.matrix(xy))
            if (!(is.null(mycall))) {
                param <- expand.grid(as.list(mycall))
                m1 <- match(names(param), names(formals(f)))
                for (i in 1:nrow(param)) {
                    formals(f)[m1] <- unclass(param[i,])
                    res[[i]] <-
                        scores.listw(nb2listw(
                            nb,
                            style = "B",
                            glist = lapply(nbdist, f), 
                            zero.policy = TRUE
                        ),
                            MEM.autocor = MEM.autocor)
                }
            }
            else {
                res[[1]] <-
                    scores.listw(nb2listw(nb, style = "B", glist = lapply(nbdist, f)),
                        MEM.autocor = MEM.autocor)
            }
        }
        else {
            res[[1]] <-
                scores.listw(nb2listw(nb, style = "B"), MEM.autocor = MEM.autocor)
        }
        res2 <-
            lapply(res, function(x)
                ortho.AIC(
                    Y = Y,
                    X = x,
                    ord.var = TRUE
                ))
        if (!(is.null(mycall))) {
            res3 <-
                data.frame(AICc = unlist(lapply(res2, function(x)
                    min(x[[1]], na.rm = TRUE))), NbVar = unlist(lapply(res2, function(x)
                        which.min(x[[1]]))))
            res3 <- cbind(param, res3)
        }
        else{
            res3 <-
                data.frame(AICc = unlist(lapply(res2, function(x)
                    min(x[[1]], na.rm = TRUE))), NbVar = unlist(lapply(res2, function(x)
                        which.min(x[[1]]))))
        }
        
        thebest <- which.min(res3$AICc)
        cat (paste("\n\nAICc for the null model:", res2[[thebest]]$AICc0, "\n"))
        cat ("\nBest spatial model:\n")
        print(res3[thebest,])
        
        return(list(all = res3, best = list(MEM = res[[thebest]], AIC = res2[[thebest]])))
        
    }
