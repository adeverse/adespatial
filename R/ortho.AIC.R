#' Compute AIC for models with orthonormal explanatory variables
#' 
#' This function compute corrected AIC for models with orthonormal and centered 
#' explanatory variables such as MEM spatial eigenfunctions. Variables are
#' sorted by their contribution to R2. .
#' 
#' It ensures that a model with k variables is the best one that can be 
#' obtained. By default, response variables are centered (model with intercept).
#' 
#' 
#' @param Y A matrix with response variables (univariate or multivariate 
#'   response)
#' @param X A set of orthonormal and centered vectors
#' @param ord.var A logical value indicating if the order of variables and 
#'   cumulative R2 must be returned
#' @return A vector with corrected AIC if \code{ord.var=FALSE}. A list if
#'   \code{ord.var=TRUE} with: \item{AICc }{Values of corrected AIC.}
#'   \item{AICc0 }{Values of corrected AIC for the null model (only intercept).}
#'   \item{ord }{Order of variables to be enter in the model} \item{R2
#'   }{Cumulative R2}
#' @author Stéphane Dray \email{stephane.dray@@univ-lyon1.fr}
#' @references Godinez-Dominguez E. and Freire J. (2003) Information-theoretic 
#'   approach for selection of spatial and temporal models of community 
#'   organization. Marine Ecology - Progress Series. 253, 17--24
#' @keywords models
#' @examples
#' 
#' y <- matrix(rnorm(50),50,1)
#' x <- svd(scale(y \%*\% c(0.1,0.5,2,0,0.7)+matrix(rnorm(250),50,5)))$u
#' res <- ortho.AIC(y,x,ord.var=TRUE)
#' minAIC <- which.min(res$AICc)
#' nvar <- length(1:minAIC)+1 # number of orthogonal vectors + 1 for intercept
#' lm1 <- lm(y~x[,res$ord[1:minAIC]])
#' summary(lm1)$r.squared # R2
#' res$R2[minAIC] # the same

#' min(res$AICc) # corrected AIC
#' extractAIC(lm1) # classical AIC
#' min(res$AICc)-2*(nvar*(nvar+1))/(nrow(x)-nvar-1) # the same
#'
#' lm2 <- lm(y~1)
#'
#' res$AICc0 # corrected AIC for the null model
#' extractAIC(lm2) # classical AIC
#' res$AICc0-2*(1*(1+1))/(nrow(x)-1-1) # the same
#'
#' @export


"ortho.AIC" <-
    function(Y, X, ord.var = FALSE) {
        # Fast Forward Selection AIC if X is orthonormal (XtX=I)
        # return a vector of AICc
        # if ord.var=TRUE, a list containing also order of variables is returned
        if (sum(apply(as.matrix(apply(X, 2, mean)), 1, function(x)
            identical(all.equal(x, 0), TRUE))) != ncol(X))
            stop("X variables are not centered")
        X <- sweep(X, 2, sqrt(colSums(X^2)), "/")
         if (!(sum(identical(all.equal(
             sum(crossprod(as.matrix(X)) - diag(ncol(X))), 0
         ), TRUE))))
             stop("X variables are not orthonormalized")
         
        f1 <- function(resp, X) {
            R2 <- t(as.matrix(X)) %*% as.matrix(Y)
            R2 <- t(R2) %*% R2
            return(sum(diag(R2)))
        }
        
        Y <- scale(Y, scale = FALSE)
        Y <- as.matrix(Y)
        R2 <- apply(X, 2, f1, resp = Y)
        SSTot <- sum(diag(t(Y) %*% Y))
        RSS <- SSTot - R2
        ordre <- order(RSS)
        RSS <- sort(RSS)
        R2 <- R2[ordre]
        RSScum <- cumsum(c(RSS[1], -R2[-1]))
        RSScum <- c(SSTot, RSScum)
        # By default, Y is centered
        # K is the number of othogonal vectors + 1 (for intercept)
        
        K <- (1 + (0:ncol(X)))
        AICtri <-
            nrow(X) * log(ifelse(RSScum <= 0, NA, RSScum / nrow(X))) + 2 * K
        correct <- 2 * (K * (K + 1)) / (nrow(X) - K - 1)
        correct <- ifelse(is.finite(correct) & (correct > 0), correct, NA)
        AICc <- AICtri + correct
        if (ord.var) {
            AICc <-
                list(
                    AICc = AICc[-1],
                    AICc0 = AICc[1],
                    ord = ordre,
                    R2 = cumsum(R2 / SSTot)
                )
        }
        if (!ord.var)
            AICc <- AICc[-1]
        return(AICc)
    }
