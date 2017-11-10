#' Function to compute eigenvectors of spatial weighting matrices
#' 
#' This function is a user-friendly way to compute eigenvectors for various 
#' possible definitions of spatial weighting matrices (W matrices). It uses the 
#' function \code{scores.listw} to build the MEM variables on the basis of the 
#' W matrix provided. Contrary to \code{\link{MEM.modsel}}, \code{test.W.R2} 
#' does not test the significance of the W matrix nor select the best subset of
#' spatial predictors within it; it simply generates and returns the complete set 
#' of MEM variables positively, negatively or both positively and negatively
#' spatially autocorrelated, depending on the specification of the argument 
#' \code{MEM.autocor}. 
#' 
#' @details This functions allows to test one binary spatial weighting matrix 
#'   (if only Y and nb are provided). It allows also to test a weighting 
#'   function based on distances (if f is provided) and a weighting function 
#'   with different values of parameters if other arguments of \code{f} are 
#'   provided. 
#'   \code{test.W.R2} is the new version of the old \code{test.W} function. The
#'   latter also generated the MEM variables of the W matrix provided (argument
#'   \code{nb}), but it additionally selected a best W matrix based on an 
#'   AIC criterion (Dray et al. 2006). This function presented two major issues:
#'   it was used without control of the global type I error rate (no correction
#'   of the significance threshold value of the W matrices), and was used in
#'   many studies to select a best subset of spatial predictors, while the method
#'   was not designed for this purpose and therefore led to highly inflated type
#'   I error rates (see Bauman et al. 2017 for details). The function 
#'   \code{test.W} should therefore not be used anymore. Instead, the user can
#'   use \code{test.W.R2} to generate the MEM variables according to different W
#'   matrices (but see warnings below), or can use the combination of 
#'   \code{\link{listw.candidates}} and \code{\link{MEM.modsel}} which is the best
#'   user-friendly way to both select a W matrix and a subset of spatial predictors
#'   within it while keeping a correct type I error rate (Bauman et al. 2018).
#'   Note that if a significance test is performed on a W matrix, the significance 
#'   threshold value must be corrected to consider the total number of W matrices 
#'   tested (Bauman et al. 2018). For instance, if the user buids W matrices using 
#'   a Gabriel's graph and a weighting function for which five parametre values are 
#'   tested, then a list of five W matrices will be returned by \code{test.W.R2}, 
#'   and the significance of each of the five W matrices should be tested at a 
#'   corrected significance threshold for five tests to avoid any type I error rate 
#'   inflation. Since \code{test.W.R2} only generates the MEM variables, the user must 
#'   test himself the significance of the global model (whole set of MEM
#'   variables generated) at the corrected significance threshold value, and
#'   then proceed to an eigenvector selection within the chosen W matrix (see 
#'   Bauman et al. 2017 for a discussion on which spatial eigenvector selection
#'   method to use in which situation).
#'   The function \code{\link{listw.candidates}} builds a list of W matrices
#'   the user is interested to compare, and this list can be entered to the
#'   function \code{\link{MEM.modsel}} which optimises the selection of the 
#'   W matrix while controling the type I error rate, and then selects the best
#'   spatial predictor subset within the optimal W matrix.
#'   
#' @param nb An object of the class \code{nb} created by functions of the 
#'   \code{spdep} package
#' @param xy Coordinates of the samples, this argument is optional and is 
#'   required only if the argument \code{f} is not null.
#' @param MEM.autocor A string indicating if all MEM must be returned or only 
#'   those corresponding to positive or negative autocorrelation; the default is 
#'   \code{positive}
#' @param style Coding scheme style (see \code{nb2listw} of the \code{spdep}
#'   package). Can take values 'W', 'B', 'C', 'U', 'minmax', and 'S'; default is
#'   'B'
#' @param f A function of the distance that can be used as a weighting spatial 
#'   function. This argument is optional
#' @param \dots Others arguments for the function \code{f}. It defines the range
#'   of parameters which will be tested
#' @return A list with the following elements: \item{param }{Only if \code{f} is
#' different than \code{NULL} and if more than one parametre value was used in
#' the function; A table of the different combinations of parametres for which a 
#' W matrix was built.} \item{MEM }{A list containing the results of 
#'   \code{scores.listw} for all the combinations of \code{param}.}
#' @author Stephane Dray and Bauman David \email{stephane.dray@@univ-lyon1.fr}, 
#' \email{davbauman@@gmail.com}
#' @seealso   \code{\link{scores.listw}}, \code{\link{listw.candidates}}, 
#' \code{\link{MEM.modsel}}
#' @references Dray S., Legendre P. and Peres-Neto P. R. (2006) Spatial 
#'   modeling: a comprehensive framework for principal coordinate analysis of 
#'   neighbor matrices (PCNM). Ecological Modelling, 196, 483--493
#'   
#' Bauman D., Drouet T., Dray S. and Vleminckx J. (2017) Disentangling good from bad 
#' practices in the selection of spatial or phylogenetic eigenvectors. Ecography 
#' *******************************
#' 
#' Bauman D., Fortin M-J, Suez M., Drouet T. and Dray S. (2018) ********************
#' 
#' @keywords spatial
#' @examples
#' 
#' if(require(ade4) & require(spdep)){
#' 
#' data(oribatid) 
#' 
#' # Generate a binary spatial weighting matrix
#' nbtri <- tri2nb(as.matrix(oribatid$xy))
#' tri.res <- test.W.R2(nbtri)
#' 
#' maxi <- max(unlist(nbdists(nbtri, as.matrix(oribatid$xy))))
#' 
#' # Generates a simple spatial weighting function of the distance
#' f1 <- function(x) {1-(x)/(maxi)}
#' tri.f1 <- test.W.R2(nbtri, f = f1, xy = as.matrix(oribatid$xy))
#' 
#' # Generates a spatial weighting function with various values of one parametre
#' f2 <- function(x, dmax, y) {1-(x^y)/(dmax)^y}
#' tri.f2 <- test.W.R2(nbtri, f = f2, y = 2:10, dmax = maxi, xy = as.matrix(oribatid$xy))
#' 
#' }
#' 
#' @importFrom spdep nb2listw
#' @export

"test.W.R2" <- function(nb, xy, style = "B", MEM.autocor = c("positive", "negative",
                                                              "all"), f = NULL, ...) 
{
  mycall <- pairlist(...)   
  res <- list()  
  MEM.autocor <- match.arg(MEM.autocor)
  control <- FALSE
  
  if (!(is.null(f))) {
    nbdist <- nbdists(nb, as.matrix(xy))
    if (!(is.null(mycall))) {   
      param <- expand.grid(as.list(mycall))
      m1 <- match(names(param), names(formals(f)))
      for (i in 1:nrow(param)) {
        formals(f)[m1] <- unclass(param[i, ])
        res[[i]] <- scores.listw(nb2listw(nb, style = style, 
                                          glist = lapply(nbdist, f), zero.policy = TRUE),
                                 MEM.autocor = MEM.autocor)
        if (i > 1) control <- TRUE
      }
    }
    else {   
      res[[1]] <- scores.listw(nb2listw(nb, style = style, glist = lapply(nbdist, f)), 
                               MEM.autocor = MEM.autocor)
    }
  }
  else {   
    res[[1]] <- scores.listw(nb2listw(nb, style = style), MEM.autocor = MEM.autocor)
  }
  
  if (control == TRUE) list(param = param, MEM = res)
  else list(MEM = res[[1]])
}