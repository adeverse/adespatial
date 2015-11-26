#' Global and local tests
#' 
#' These two Monte Carlo tests are used to assess the existence of 'global' and
#' 'local' spatial structures, corresponding respectively to positive and
#' negative Moran's I .\cr
#' 
#' They rely on the decomposition of a data matrix X into global and local
#' components using multiple regression on Moran's Eigenvector Maps (MEMs).
#' They require a data matrix (X) and a list of weights derived from a
#' connection network. X is regressed onto global MEMs (U+) in the global test
#' and on local ones (U-) in the local test. One mean \eqn{R^2}{R^2} is
#' obtained for each MEM, the k highest being summed to form the test
#' statistic.
#' 
#' The reference distribution of these statistics are obtained by randomly
#' permuting the rows of X.
#' 
#' These tests were originally part of the adegenet package for R.
#' 
#' @aliases global.rtest local.rtest
#' @param X a data matrix, with variables in columns
#' @param listw a list of weights of class \code{listw}. Can be obtained easily
#' using the function \code{chooseCN}.
#' @param k integer: the number of highest \eqn{R^2}{R^2} summed to form the
#' test statistics
#' @param nperm integer: the number of randomisations to be performed.
#' @return An object of class \code{randtest}.
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @references Jombart, T., Devillard, S., Dufour, A.-B. and Pontier, D.
#' Revealing cryptic spatial patterns in genetic variability by a new
#' multivariate method. \emph{Heredity} doi: 10.1038/hdy.2008.34.
#' @keywords multivariate spatial
#' @examples
#' 
#' 
#' # wait for a generic dataset
#' 
#' 
#' @importFrom ade4 as.randtest
#' @export global.rtest local.rtest

global.rtest <- function(X, listw, k = 1, nperm = 499){
    if (!inherits(listw, "listw")) stop("object of class 'listw' expected")
    if (listw$style != "W") stop("object of class 'listw' with style 'W' expected")
    if(any(is.na(X))) stop("NA entries in X")
    
    n <- nrow(X)
    X <- scalewt(X)
    
    # computation of U+
    temp <- orthobasis.listw(listw)
    val <- attr(temp,"values")
    U <- as.matrix(temp)
    Upos <-  U[,val > -1/(n-1)]
    
    # test statistic
    calcstat <- function(X, k){
        R <- (t(X) %*% Upos) / n
        R2 <- R*R
        temp <- sort(apply(R2, 2, mean), decreasing=TRUE)
        stat <- sum(temp[1:k])
        return(stat)
    }
    
    ini <- calcstat(X, k)
    
    sim <- sapply(1:nperm, function(i) calcstat(X[sample(1:n),], k ))
    
    res <- as.randtest(sim = sim, obs = ini, alter = "greater")
    res$call <- match.call()
    
    return(res)
} #end global.rtest

local.rtest <- function(X, listw, k = 1, nperm = 499){
    if (!inherits(listw, "listw")) stop("object of class 'listw' expected")
    if (listw$style != "W") stop("object of class 'listw' with style 'W' expected")
    if(any(is.na(X))) stop("NA entries in X")
    
    n <- nrow(X)
    X <- scalewt(X)
    
    # computation of U-
    temp <- orthobasis.listw(listw)
    val <- attr(temp,"values")
    U <- as.matrix(temp)
    Uneg <-  U[,val < -1/(n-1)]
    
    X <- scalewt(X)
    
    # test statistic
    calcstat <- function(X, k){
        R <- (t(X) %*% Uneg) / n
        R2 <- R*R
        temp <- sort(apply(R2, 2, mean),decreasing=TRUE)
        stat <- sum(temp[1:k])
        return(stat)
    }
    
    ini <- calcstat(X, k)
    
    sim <- sapply(1:nperm, function(i) calcstat(X[sample(1:n),], k ))
    
    res <- as.randtest(sim = sim, obs = ini, alter = "greater")
    res$call <- match.call()
    
    return(res)
} 

