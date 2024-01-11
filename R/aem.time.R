#' AEM for time series
#'
#' This function constructs AEM eigenfunctions for multi-scale analysis of a
#' regular time series or spatial transect of univariate or multivariate data.
#'
#' @param n Numeric. Number of points in the series.
#' @param w A vector of weights to be applied to the edges (columns of matrix
#'   E). Equal weights are used if no vector \code{w} is provided. The length of
#'   vector \code{w} must be (\code{n}-1) where \code{n} is the number of points
#'   in the spatial or temporal series.
#' @param moran Logical. If \code{TRUE}, Moran's I are computed for all AEM. If
#'   \code{FALSE} (default value), Moran's I are not computed.
#'
#' @details
#'
#' Time series represent a form of directional stochastic process. To emphasize
#' the directional nature of the process influencing the data, AEM analysis,
#' which was designed to take trends into account, should be applied to the
#' non-detrended series. MEM analysis (see \code{scores.listw}) can be applied
#' to data series that were detrended to remove the directional component as
#' recommended by Blanchet et al. (2008, 2011) and  Legendre & Legendre (2012,
#' Subsection 14.1.2). Detrended palaeoecological sediment core data, for
#' example, could be studied by MEM analysis.
#'
#' No data file needs to be provided to this function. The AEM eigenvectors are
#' constructed from a matrix E generated from the regular sequence of points
#' along the series.
#'
#' A vector of weights \code{w} can be provided, representing the ease of
#' communication of matter, energy or information among the points. The most
#' simple form would be the inverse of (d/dmax) where d is the distance between
#' adjacent nodes and dmax is the maximum distance between adjacent nodes in the
#' spatial or time series. More general forms of weights may represent the
#' inverse of landscape resistance to the movement of organisms, propagules,
#' genes, etc.
#'
#' If the calculation of Moran's I is requested, the point coordinates are
#' generated from the point positions along the series.
#'
#' @return
#'
#' \item{E}{Nodes-by-edges matrix E. }
#' \item{values}{Eigenvalues of the
#' principal component analysis of E. }
#' \item{aem}{Matrix of AEM eigenfunctions
#' normalized to unit length. }
#' \item{Moran}{Moran's I statistics tested by a bilateral test with 999 permutations}
#' \item{listw}{An object of class \code{listw} with the associated spatial weighting matrix}
#'
#' @author Pierre Legendre and F. Guillaume Blanchet
#'
#' @seealso \code{\link{aem}}, \code{scores.listw}
#'
#' @references
#'
#' Blanchet F.G., P. Legendre and Borcard D. (2008) Modelling directional
#' spatial processes in ecological data. \emph{Ecological Modelling}, 215,
#' 325-336.
#'
#' Blanchet F.G., P. Legendre, R. Maranger, D. Monti, and P. Pepin. (2011)
#' Modelling the effect of directional spatial ecological processes at different
#' scales. \emph{Oecologia}, 166, 357-368.
#'
#' Legendre, P. and L. Legendre (2012) \emph{Numerical Ecology}, 3rd English
#' edition. Elsevier Science BV, Amsterdam.
#'
#' Legendre, P. and O. Gauthier (2014) Statistical methods for temporal and
#' space-time analysis of community composition data. \emph{Proceedings of the
#' Royal Society B - Biological Sciences}, 281, 20132728.
#'
#' @examples
#'
#' # Time series containing 20 equispaced observations
#' out <- aem.time(20, moran = TRUE)
#'
#' # Time series containing 20 observations with unequal spacing
#' # Generate (n-1) random interpoint distances
#' distances <- runif(19,1,5)
#'
#' # Compute weights representing the ease of communication among points
#' w <- 1/(distances/max(distances))
#'
#' # Compute the AEM eigenfunctions
#' out <- aem.time(20, w = w, moran = TRUE)
#'
#' @keywords multivariate
#' @keywords spatial
#' @export

aem.time <- function(n,
    w = NULL,
    moran = FALSE) {
    #
    # Construct AEM eigenfunctions for a regular time series
    # n = number of points
    # w = vector of weights. The weights can be the inverse of the interval lengths if the observations are not regularly spaced. Default: equal weights.
    # moran: compute Moran's I for each AEM eigenfunction
    #
    # Authors: Pierre Legendre and F. Guillaume Blanchet, March 2012
    
    # Normalize a vector (to length 1)
    normalize <- function(vec)
        vec / sqrt(sum(vec ^ 2))
    ###  End internal functions
    #
    epsilon <- sqrt(.Machine$double.eps)
    #
    # Construct matrix E
    E <- matrix(0, n, (n - 1))
    rownames(E) <- paste("site", 1:n, sep = ".")
    colnames(E) <- paste("E", 1:(n - 1), sep = "")
    for (i in 2:n)
        E[i, 1:(i - 1)] <- 1
    #
    # Apply weights if provided
    if (!is.null(w)) {
        if (length(w) != (i - 1))
            stop("Length of vector w not equal to (n-1)")
        E <- E %*% diag(w)
    } else {
        w <- rep(1, n - 1)
    }
    #
    # Compute AEM eigenfunctions
    E.c <- scale(E, center = TRUE, scale = FALSE)
    E.svd <- svd(E.c)
    k <- length(which(E.svd$d > epsilon))
    # Normalize the AEM eigenfunctions
    E.svd$u[, 1:k] <- apply(E.svd$u[, 1:k], 2, normalize)
    
    mat01 <- matrix(0, n, n)
    for (i in 1:(n - 1))
        mat01[i, i + 1] <- w[i]
    lw <- mat2listw(mat01, style = "B", zero.policy = TRUE)
    
    out <- list(E = E,
        values = E.svd$d[1:k] ^ 2 / (n - 1),
        aem = E.svd$u[, 1:k], listw = lw)
    
    if (moran) 
        out$Moran <- moran.randtest(E.svd$u[,1:k], alter = "two-sided", listw = lw)
  
    out
}
