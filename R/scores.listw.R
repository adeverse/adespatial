#' Function to compute and manage Moran's Eigenvector Maps (MEM) of a listw
#' object
#'
#' These functions compute MEM (i.e., eigenvectors of a doubly centered spatial
#' weighting matrix). Corresponding eigenvalues are linearly related to Moran's
#' index of spatial autocorrelation.
#'
#' Testing the nullity of eigenvalues is based on E(i)/E(1) where E(i) is i-th
#' eigenvalue and E(1) is the maximum absolute value of eigenvalues
#'
#' @aliases scores.listw mem orthobasis.listw [.orthobasisSp
#' @param listw An object of the class \code{listw} created by functions of the
#'   \code{spdep} package
#' @param wt A vector of weights. It is used to orthogonalize the eigenvectors.
#'   It could be useful if MEM are used in weighted regression or canonical
#'   correspondence analysis
#' @param MEM.autocor A string indicating if all MEMs must be returned or only
#'   those corresponding to non-null, positive or negative autocorrelation. The
#'   difference between options \code{all} and \code{non-null} is the following:
#'   when there are several null eigenvalues, option \code{all} removes only one
#'   of the eigenvectors with null eigenvalues and returns (n-1) eigenvectors,
#'   whereas \code{non-null} does not return any of the eigenvectors with null
#'   eigenvalues.
#' @param store.listw A logical indicating if the spatial weighting matrix
#'   should be stored in the attribute \code{listw} of the returned object
#' @param x An object of class \code{orthobasisSp}.
#' @param i,j Elements to extract (integer or empty): index of rows (i) and
#'   columns (j).
#' @return An object of class \code{orthobasisSp} , subclass \code{orthobasis}.
#'   The MEMs are stored as a \code{data.frame}. It contains several attributes
#'   (see \code{?attributes}) including: \itemize{\item \code{values}: The
#'   associated eigenvalues. \item \code{listw}: The associated spatial
#'   weighting matrix (if \code{store.listw = TRUE}). }
#' @author Stéphane Dray \email{stephane.dray@@univ-lyon1.fr}
#' @seealso \code{\link[spdep]{nb2listw}} \code{\link[ade4]{orthobasis}}
#' @references Dray, S., Legendre, P., and Peres-Neto, P. R. (2006). Spatial
#'   modeling: a comprehensive framework for principal coordinate analysis of
#'   neighbor matrices (PCNM). \emph{Ecological Modelling} \bold{196}, 483--493.
#'
#'   Griffith D. A. (1996) Spatial autocorrelation and eigenfunctions of the
#'   geographic weights matrix accompanying geo-referenced data. \emph{Canadian
#'   Geographer} \bold{40}, 351--367.
#' @keywords spatial
#' @examples
#'
#' if(require("ade4", quietly = TRUE) & require("spdep", quietly = TRUE)){
#' data(oribatid)
#' nbtri <- tri2nb(as.matrix(oribatid$xy))
#' sc.tri <- scores.listw(nb2listw(nbtri, style = "B"))
#' summary(sc.tri)
#' }
#' if(require("adegraphics", quietly = TRUE)){
#' s.value(oribatid$xy,sc.tri[,1:9])
#' plot(sc.tri[,1:6], oribatid$xy, pSp.cex = 5, pSp.alpha = 0.5, pbackground.col = 'lightblue')
#' }
#'
#' @export scores.listw orthobasis.listw mem
#' @importFrom ade4 bicenter.wt
#' @importFrom spdep listw2mat
#' @rdname mem

scores.listw <- function (listw,  wt = rep(1, length(listw$neighbours)), MEM.autocor = c("non-null", "all", "positive", 
                                                                                         "negative"), store.listw = FALSE) 
{
    if (!inherits(listw, "listw")) 
        stop("not a listw object")
    MEM.autocor <- match.arg(MEM.autocor)
    
    wt <- wt / sum(wt)
    w <- listw2mat(listw)
    sumW <- sum(w)
    n <- nrow(w)
    symmetric <- isSymmetric.matrix(w, check.attributes = FALSE)
    if (symmetric == FALSE) {
        w <- (w + t(w))/2
    }
    w <- bicenter.wt(w, row.wt = wt, col.wt = wt)
    wtsqrt <- sqrt(wt)
    w <- w * wtsqrt
    w <- t(t(w) * wtsqrt)
    res <- eigen(w, symmetric = TRUE)
    eq0 <- which(apply(as.matrix(res$values/max(abs(res$values))), 
                       1, function(x) identical(all.equal(x, 0), TRUE)))
    if (length(eq0) == 0) {
        stop("Illegal matrix: no null eigenvalue")
    }
    
    if(MEM.autocor == "all"){
        if (length(eq0) == 1) {
            res$values <- res$values[-eq0]
            res$vectors <- res$vectors[, -eq0]
            
        } else if (length(eq0) > 1) {
            w <- cbind(wt, res$vectors[, eq0])
            w <- qr.Q(qr(w))
            res$values[eq0] <- 0
            res$vectors[, eq0] <- w[, -ncol(w)]
            res$values <- res$values[-eq0[1]]
            res$vectors <- res$vectors[, -eq0[1]]
        }
    } else if(MEM.autocor == "non-null"){        
        res$values <- res$values[-eq0]
        res$vectors <- res$vectors[, -eq0]
    } else if (MEM.autocor == "positive") {
        posi <- 1:(min(eq0)-1)
        res$values <- res$values[posi]
        res$vectors <- res$vectors[, posi]
    } else if (MEM.autocor == "negative") {
        neg <- (max(eq0) + 1):ncol(res$vectors)
        res$values <- res$values[neg]
        res$vectors <- res$vectors[, neg]
    }
    
    a0 <- as.data.frame(res$vectors) / wtsqrt
    z <- res$values
    row.names(a0) <- attr(listw,"region.id")
    names(a0) <- paste("MEM", 1:ncol(a0), sep = "")
    attr(a0,"values") <- z
    attr(a0,"weights") <- wt
    if(store.listw)
        attr(a0,"listw") <- listw
    attr(a0,"call") <- match.call()
    attr(a0,"class") <- c("orthobasisSp", "orthobasis","data.frame")
    
    return(a0)
}

#' @rdname mem
mem <- function (listw,  wt = rep(1, length(listw$neighbours)), 
                 MEM.autocor = c("non-null", "all", "positive", "negative"), store.listw = FALSE) {
    res <- scores.listw(listw = listw, wt = wt, MEM.autocor = MEM.autocor, store.listw = store.listw)
    attr(res,"call") <- match.call()
    return(res)
}

#' @rdname mem
orthobasis.listw <- function (listw,  wt = rep(1, length(listw$neighbours)), 
                              MEM.autocor = c("non-null", "all", "positive", "negative"), store.listw = FALSE) {
    res <- scores.listw(listw = listw, wt = wt, MEM.autocor = MEM.autocor, store.listw = store.listw)
    attr(res,"call") <- match.call()
    return(res)
}

#' @rdname mem
#' @export
"[.orthobasisSp" <- function (x, i, j) {
    
    ## i: index of rows
    ## j: index of columns
    res <- as.data.frame(x)
    if(!missing(i)){
        res <- res[i, , drop = FALSE]
        attr(res,"values") <- attr(x,"values")
        attr(res,"weights") <- attr(x,"weights")[i]
    }
    if(!missing(j)){
        res <- res[, j, drop = FALSE]
        attr(res,"values") <- attr(x,"values")[j]
        attr(res,"weights") <- attr(x,"weights")
    }
    attr(res,"call") <- match.call()
    attr(res,"class") <- attr(x,"class") 
    if(ncol(res) == 1)
        res <- as.vector(as.matrix(res))
    return(res)
}
