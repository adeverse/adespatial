#' Moran spectral randomization
#' 
#' This function allows to generate spatially-constrained random variables
#' preserving the global autocorrelation (Moran's I) and the spatial structures
#' at multiple scales. Multiscale property is defined by the power spectrum
#' (i.e. decomposition of the variance of the original variables) on a basis of
#' orthonormal eigenvectors (Moran's Eigenvector Maps, MEM). The function
#' provides methods for univariate randomization, joint randomization of a
#' group of variables while keeping within-group correlations fixed and
#' univariate randomization with a fixed correlation between original data and
#' randomized replicates.
#' 
#' Three procedures are implemented in the function. The "pair" procedure is
#' the more general as it can be applied in the three cases (univariate,
#' univariate with fixed correlation and multivariate). This procedure
#' preserves the power spectrum by pair of MEMs but not strictly the global
#' autocorrelation level (Moran's I). The "singleton" procedure can be used for
#' univariate and multivariate cases. It preserves strictly the global level of
#' autocorrelation and the power spectrum. The "triplet" procedure can only be
#' applied in the univariate case. It preserves the power spectrum by triplet
#' of MEMs and strictly the global autocorrelation level.
#' 
#' @param x a \code{vector}, a \code{matrix} or a \code{data.frame} with the
#' original variables. If \code{NCOL(x) > 1}, then the joint randomization
#' procedure that preserves the correlations among variables is used.
#' @param listwORorthobasis an object of the class \code{listw} (spatial
#' weights) created by the functions of the \pkg{spdep} package or an object of
#' class \code{orthobasis}
#' @param nrepet an \code{integer} indicating the number of replicates
#' @param method an character specifying which algorithm should be used to
#' produce spatial replicates (see Details).
#' @param cor.fixed if not missing, the level of correlation between the
#' original variable and its randomized replicates
#' @param nmax the number of trials used in the "triplet" procedure.
#' @param simplify A logical value. If \code{TRUE}, the outputs for univariate
#' procedures are returned in a matrix where each column corresponds to a
#' replicate. If \code{FALSE}n a \code{list} is returned.
#' @return Either a matrix (if \code{simplify} is \code{TRUE}) or a list with
#' randomized replicates.
#' @author Stephane Dray \email{stephane.dray@@univ-lyon1.fr} and Helene H
#' Wagner \email{helene.wagner@@utoronto.ca}
#' @seealso \code{\link{scores.listw}}, \code{\link[spdep]{nb2listw}}
#' @references Wagner, H.H. and Dray S. (2015) Generating spatially-constrained
#' null models for irregularly spaced data using Moran spectral randomization
#' methods
#' @keywords spatial
#' @examples
#' 
#' library(spdep)
#' x1 <- matrix(rnorm(81*5), nrow = 81)
#' lw1 <- nb2listw(cell2nb(9, 9))
#' 
#' moran.mc(x1[,1], lw1, 2)$statistic
#' 
#' ## singleton
#' x1.1 <- msr(x1[,1], lw1, nrepet = 9, method = "singleton")
#' apply(x1.1, 2, function(x) moran.mc(x, listw = lw1, nsim = 2)$statistic)
#' 
#' ## triplet
#' x1.2 <- msr(x1[,1], lw1, nrepet = 9, method = "triplet")
#' apply(x1.2, 2, function(x) moran.mc(x, listw = lw1, nsim = 2)$statistic)
#' 
#' ## pair
#' x1.3 <- msr(x1[,1], lw1, nrepet = 9, method = "pair")
#' apply(x1.3, 2, function(x) moran.mc(x, listw = lw1, nsim = 2)$statistic)
#' 
#' ## pair with cor.fixed
#' x1.4 <- msr(x1[,1], lw1, nrepet = 9, cor.fixed = 0.5)
#' apply(x1.4, 2, function(x) moran.mc(x, listw = lw1, nsim = 2)$statistic)
#' cor(x1[,1], x1.4)
#' 
#' ## pair preserving correlations for multivariate data
#' x1.5 <- msr(x1, lw1, nrepet = 9, cor.fixed = 0.5)
#' cor(x1)
#' lapply(x1.5, cor)
#' 
#' apply(x1, 2, function(x) moran.mc(x, listw = lw1, nsim = 2)$statistic)
#' apply(x1.5[[1]], 2, function(x) moran.mc(x, listw = lw1, nsim = 2)$statistic)
#' 
#' ## singleton preserving correlations for multivariate data
#' x1.6 <- msr(x1, lw1, nrepet = 9, method = "singleton")
#' cor(x1)
#' lapply(x1.6, cor)
#' 
#' apply(x1, 2, function(x) moran.mc(x, listw = lw1, nsim = 2)$statistic)
#' apply(x1.6[[1]], 2, function(x) moran.mc(x, listw = lw1, nsim = 2)$statistic)
#' 
#' @export msr
msr <- function(x, listwORorthobasis, nrepet = 99, method = c("pair", "triplet", "singleton"), cor.fixed, nmax = 100, simplify = TRUE){

    if (inherits(listwORorthobasis, "listw")){
        tmp <- scores.listw(listwORorthobasis, MEM.autocor = "all")
        mem <- list(vectors = as.matrix(tmp), values = attr(tmp, "values"))
        
    } else if(inherits(listwORorthobasis, "orthobasis")){
          mem <- list(vectors = as.matrix(listwORorthobasis), values = attr(listwORorthobasis, "values"))
          if(!identical(attr(listwORorthobasis, "weights"), rep(1/nrow(mem$vectors), nrow(mem$vectors))))
              stop("Only implemented for uniform weights")
          if(ncol(mem$vectors) != NROW(x) - 1)
              stop("orthobasis object should contain n-1 vectors")
      } else {
            stop("listwORorthobasis is not a listw or orthobasis object")
        }

    mem$vectors <-  mem$vectors / sqrt(nrow(mem$vectors))
    r <- cor(as.matrix(x), mem$vectors) 
    x.mean <- colMeans(as.matrix(x))
    x.sd <- apply(as.matrix(x), 2, sd)
    method <- match.arg(method)
    res <- vector("list", nrepet) ## results are stored in a list where each elemet corresponds to a repetion
    
    if(!missing(cor.fixed)){
        method <- "pair"
    }
    
     if((NCOL(x) > 1) & method == "triplet"){
         stop("The 'triplet' method can not handle multivariate data. Use 'pair' or 'singleton' instead")
     }
    
    for(nr in 1:nrepet){
        ## singleton
        if(method == "singleton"){
          
          if(NROW(r) == 1)
              a <- sapply(r, genbysingleton)
          else
              a <- t(apply(r, 2, genbysingleton)) 
          xnew <- mem$vectors%*%a
        }

        ## pair
        if(method == "pair"){
            if(missing(cor.fixed))
                cor.fixed <- NULL
            pairs <- cutbypair(mem$values)
            a <- array(NA, c(nrow(pairs), 2, NCOL(x)))
            for(i in 1:nrow(pairs)){
                if(sum(is.na(pairs[i,])) == 0){
                    a[i,,] <- genbypair(r[, pairs[i,], drop = FALSE], cor.fixed = cor.fixed) 
                } else {
                    a[i,1,] <- genbysingleton(r[, pairs[i,1], drop = FALSE])
                }
            }
            xnew <- mem$vectors[, t(pairs)[1:length(mem$values)]] %*% apply(a, 3, function(x) t(x)[1:length(mem$values)])
        }
        
        ## triplet
        if(method == "triplet"){
            triplets <- cutbytriplet(mem$values, method = "within")
            a <- matrix(NA, nrow = nrow(triplets), ncol = 3)
            generateby <- rep("T", length(r))
            
            for(i in 1:nrow(triplets)){
                
                if(sum(is.na(triplets[i,])) == 0){
                    niter <- 1
                    repeat{
                        a[i,] <- genbytriplet(mem$values[triplets[i,]], r[triplets[i,]])
                        if(sum(is.na(a[i,])) == 0) ## ok
                            niter <- nmax + 1
                        else
                            niter <-  niter + 1 ## try again
                        if(niter >= nmax)
                            break
                    }
                    
                } else if(sum(is.na(triplets[i,])) == 1){
                    a[i,1] <- genbysingleton(r[triplets[i, 1]])
                    a[i,2] <- genbysingleton(r[triplets[i, 2]])       
                 } else if(sum(is.na(triplets[i,])) == 2)
                      a[i,1] <- genbysingleton(r[triplets[i, 1]])
                
            }
            
            ## fill the NA with the singleton method (could be by pair)
            idxNA <- triplets[(is.na(a))]
            idxNA <- idxNA[!is.na(idxNA)]
            
            if(length(idxNA) > 0){
                a[is.element(triplets,idxNA)] <- sapply(idxNA, function(i) genbysingleton(r[i]))
                generateby[idxNA] <- "S"
            }
            
            ## generate the new variable
            xnew <- mem$vectors[, t(triplets)[1:length(mem$values)]]%*%t(a)[1:length(mem$values)]
            attr(xnew, "method") <- as.factor(generateby)           
        }
        
        xnew <- sweep(sweep(sqrt(NROW(x)-1)*xnew, 2, x.sd, "*"), 2, x.mean, "+")
        res[[nr]] <- xnew
    }

    ## for univariate x, returns a matrix instead of a list
    if((NCOL(x) == 1) & simplify){
        res <- do.call("cbind", res)
    }
    
    return(res)
}


genbysingleton <- function(r){
    ## r is the matrix (nvar x 1) of correlation with MEM
    ## The functions returns a new value of r (stored in a) so that the contributions to R2 and I are preserved
    ## this is simply performed by switching the signs in the same manner for all variables
    return(sample(c(-1, 1), 1) *  r)
}

genbypair <- function(r, cor.fixed = NULL){
    ## r is a matrix (nvar by 2) with correlation of each variable with a pair of MEMs
    ## cor.fixed is (if not NULL) the value of the correlation with the original variable to be preserved (only for univariate case)
    ## The functions returns new values of r (stored in a) so that the contributions to R2 is preserved (but not I)

    nvar <- nrow(r)
    a <- matrix(0, nrow = nvar, ncol = 2)
    R2 <- rowSums(r^2) ## contribution to global R2 
    
    ##  phi: draws angle from uniform distribution
    if(nvar > 1){
        phi <- atan2(r[,2], r[,1]) + runif(1, 0, 2*pi)
    } else {
        if(is.null(cor.fixed)){
            phi <- runif(1, -pi, pi)
        } else {
            phi <- atan2(r[,2], r[,1]) + sample(c(-1, 1), 1) * acos(cor.fixed)
        }
    }

    ## determine the new values for r (i.e. a)
    R <- sqrt(R2)
    a[,1] <- R * cos(phi)
    a[,2] <- R * sin(phi) 
    
    return(as.vector(t(a)))
}



genbytriplet <- function(lambda, r){
    ## lambda is a vector (3 by 1) with the eigenvalues associated to MEMs
    ## r is a vector (3 by 1) giving the correlation with each MEM
    ## The functions returns new values of r (stored in a) so that the contributions to R2 and I are preserved
      
    a <- rep(NA, 3)
    r2 <- r^2
    R2 <- sum(r2) ## contribution to global R2 
    I <- sum(r2 * lambda) ## contribution to global I
     
    ##  theta: draws first angle from uniform distribution
    theta <- runif(1, -pi, pi)
   
    ## phi: find the second angle preserving the contributions
    sin.phi.squared <- (I/R2 - lambda[3] - sin(theta)^2 * 
                          (lambda[1] - lambda[3])) / (sin(theta)^2 * (lambda[2] - lambda[1]))

    if(is.finite(sin.phi.squared)){
        if((sin.phi.squared <= 1) & (sin.phi.squared > 0)){
            ## existing condition
            phi <- asin(sqrt(sin.phi.squared))
            rndsign <- function() sample(c(-1, 1), 1) 
            phi <- rndsign() * phi
            
            ## determine the new values for r (i.e. a)
            R <- sqrt(R2)
            a[1] <- rndsign() * R * cos(phi) * sin(theta)
            a[2] <- rndsign() * R * sin(phi) * sin(theta)
            a[3] <- rndsign() * R * cos(theta)
        }
    }
    
    return(a)
}




cutbypair <- function(lambda, method = c("consecutive")){
    ## determine the indexes for pairs using lambda values
    ## only 'consecutive' is implemented
    
    n <- length(lambda)
    n.pair <- n%/%2
    mod <- n%%2
    perm <- 1:n
    idx0 <- sample(perm, 1)
    if(mod > 0)
        perm <- perm[-idx0]
    res <- matrix(perm, byrow = TRUE, nrow = n.pair)
    if(mod > 0){
        res <- rbind(res, c(idx0, NA))
    }
    
    return(res)
}


cutbytriplet <- function(lambda, method = c("consecutive", "random", "within")){
    ## determine the indexes for triplet using lambda values
    ## in method 'consecutive', we force that non-triplet corresponds to smaller eigenvalues to minimize biases
    method <- match.arg(method)
    n <- length(lambda)
    n.triplet <- n%/%3
    mod <- n%%3
    if(method == "within"){
        npos <- sum(lambda > 0)
        nneg <- n - npos
        mod.pos <- npos%%3
        mod.neg <- nneg%%3
        perm.pos <- sample(npos) + nneg
        perm.neg <- sample(nneg)
        if(mod == 0){ ## three cases (0,0), (1,2) or (2,1)
            perm <- c(perm.pos, perm.neg)
        } else if(mod == 1){ ## three cases (1,0), (0,1) or (2,2)
            if(mod.pos == 1){
                perm <- c(perm.pos[-1], perm.neg, perm.pos[1])
            } else {
                perm <- perm <- c(perm.pos, perm.neg)
            } 
        } else if(mod == 2){ ## three cases (1,1), (2,0) or (0,2)
            if(mod.pos == 1){
                perm <- c(perm.pos[-1], perm.neg, perm.pos[1])
            } else if(mod.pos == 2){
                perm <- c(perm.pos[-c(1,2)], perm.neg, perm.pos[c(1,2)])
            } else if(mod.pos == 0){
                perm <- c(perm.pos, perm.neg)
            } 
        }
            
        res <- matrix(perm[1:(n.triplet * 3)], byrow = TRUE, nrow = n.triplet)
        if(mod > 0){
            res <- rbind(res, c(perm[((n.triplet * 3) + 1) : length(perm)], rep(NA, 3 - mod)))
        }
    } else if(method == "random"){
        perm <- sample(n)
        res <- matrix(perm[1:(n.triplet * 3)], byrow = TRUE, nrow = n.triplet)
        if(mod > 0){
            res <- rbind(res, c(perm[((n.triplet * 3) + 1) : length(perm)], rep(NA, 3 - mod)))
        }
    } else if(method == "consecutive"){
        perm <- 1:n
        idx0 <- which(rank(abs(lambda), ties.method = "random") <= mod)
        if(mod > 0)
            perm <- perm[-idx0]
        res <- matrix(perm, byrow = TRUE, nrow = n.triplet)
        if(mod > 0){
            res <- rbind(res, c(idx0, rep(NA, 3 - mod)))
        }
    }
    return(res)
}





