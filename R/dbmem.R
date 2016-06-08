#'dbMEM spatial eigenfunctions
#'
#'Compute distance-based Moran's eigenvector maps (dbMEM, also called dbMEM 
#'spatial eigenfunctions) from a geographic distance matrix, in view of spatial 
#'eigenfunction analysis.
#'
#'@param xyORdist Either a matrix of spatial coordinates or a distance matrix 
#'  (class \code{dist}).
#'@param thresh A threshold value for truncation of the geographic distance 
#'  matrix. If \code{thresh=NULL}, the length of the longest edge of the minimum
#'  spanning tree will be used as the threshold (as returned by the function 
#'  \code{give.thresh}).
#'@param MEM.autocor A string indicating if all MEMs must be returned or only 
#'  those corresponding to non-null, positive or negative autocorrelation
#'@param store.listw A logical indicating if the spatial weighting matrix should
#'  be stored in the attribute \code{listw} of the returned object
#'@param silent A logical indicating if some information should be printed 
#'  during computation
#'  
#'@return An object of class \code{orthobasisSp} , subclass \code{orthobasis}. 
#'  The dbMEM eigenfunctions (principal coordinates of the truncated distance 
#'  matrix) are stored as a \code{data.frame}. It contains several attributes 
#'  (see \code{?attributes}) including: \itemize{\item \code{values}: The dbMEM 
#'  eigenvalues. \item \code{listw}: The associated spatial weighting matrix (if
#'  \code{store.listw = TRUE}). }
#'  
#'@details dbMEM eigenfunctions were called PCNM in early papers (Borcard and 
#'  Legendre 2002, Borcard et al. 2004). There is a small difference in the 
#'  computation: to construct PCNMs, the distance matrix subjected to PCoA 
#'  contained zeros on the diagonal. In dbMEM, the matrix contains 4*thresh 
#'  values on the diagonal. The result is that the dbMEM eigenvalues are smaller
#'  than the PCNM eigenvalues by a constant (equal to (n.sites *
#'  (4*thresh)^2)/2). The dbMEM eigenvalues are proportional to Moran's I
#'  coefficient of spatial correlation (Dray et al. 2006; Legendre and Legendre
#'  2012). The dbMEM eigenvectors only differ from the PCNM eigenvectors by a
#'  multiplicative constant; this has no impact on the use of MEMs as
#'  explanatory variables in linear models. In this implementation, dbMEM
#'  eigenvectors have a norm equal to 1 (using the uniform weigts 1/n.sites).
#'  
#'  
#'  If a truncation value is not provided, the largest distance in a minimum 
#'  spanning tree linking all sites on the map is computed (returned by the 
#'  function \code{give.thresh}). That value is used as the truncation threshold
#'  value (thresh).
#'  
#'  
#'@author  Stéphane Dray \email{stephane.dray@@univ-lyon1.fr}, Pierre Legendre, 
#'  Daniel Borcard and F. Guillaume Blanchet
#'  
#'@references
#'
#'Borcard, D. and P. Legendre. 2002. All-scale spatial analysis of ecological 
#'data by means of principal coordinates of neighbour matrices. Ecological 
#'Modelling 153: 51-68.
#'
#'Borcard, D., P. Legendre, C. Avois-Jacquet and H. Tuomisto. 2004. Dissecting 
#'the spatial structure of ecological data at multiple scales. Ecology 85: 
#'1826-1832.
#'
#'Dray, S., P. Legendre and P. R. Peres-Neto. 2006. Spatial modelling: a 
#'comprehensive framework for principal coordinate analysis of neighbour 
#'matrices (PCNM). Ecological Modelling 196: 483-493.
#'
#'Legendre, P. and L. Legendre. 2012. Numerical ecology, 3rd English edition. 
#'Elsevier Science BV, Amsterdam.
#'
#'@seealso \code{\link{give.thresh}}, \code{\link{mem}}
#'  
#' @examples
#' if(require("ade4", quietly = TRUE) & require("adegraphics", quietly = TRUE)){
#'
#' data(oribatid)
#' mite <- oribatid$fau      # 70 peat cores, 35 species
#' mite.xy <- oribatid$xy    # Geographic coordinates of the 70 cores
#'
#'
#' # thresh=1.012 is the value used in Borcard and Legendre (2002)
#' mite.dbmem1 <- dbmem(mite.xy, thresh=1.012)
#' mite.dbmem1
#' 
#' # Plot the associated spatial weighting matrix
#' s.label(mite.xy, nb = attr(mite.dbmem1, "listw"))
#'
#' # Plot maps of the first 3 dbMEM eigenfunctions
#' s.value(mite.xy, mite.dbmem1[,1:3])
#'
#' # Compute and test associated Moran's I values
#' # Eigenvalues are proportional to Moran's I 
#'
#' test <- moran.randtest(mite.dbmem1, nrepet = 99)
#' plot(test$obs, attr(mite.dbmem1, "values"), xlab = "Moran's I", ylab = "Eigenvalues")
#' 
#' # Decreasing values of Moran's I for the successive MEM. 
#' # The red line is the expected value of Moran's I under H0.
#' 
#' plot(test$obs, xlab="MEM rank", ylab="Moran's I")
#' abline(h=-1/(nrow(mite.xy) - 1), col="red")
#'
#' # Compute only the dbmem with positive eigenvalues (and positive Moran's I)
#' mite.dbmem2 <- dbmem(mite.xy, thresh=1.012, MEM.autocor="positive")
#' # or:  mite.dbmem2 <- dbmem(dist(mite.xy), thresh=1.012, MEM.autocor="positive")
#' mite.dbmem2
#' 
#' # Examine the eigenvalues
#' attributes(mite.dbmem2)$values 
#' # or:  attr(mite.dbmem2, "values")
#' 
#' # Examine (any portion of) the dbmem spatial eigenvectors
#' tmp <- as.matrix(mite.dbmem2)
#' tmp[1:10,1:6]
#'}
#'
#'@importFrom ade4 dudi.pco
#'@importFrom stats dist
#'@importFrom spdep dnearneigh nbdists nb2listw
#'@export dbmem
#'  
#'  

'dbmem' <-
    function(xyORdist,
             thresh = NULL,
             MEM.autocor = c("non-null", "all", "positive",
                             "negative"),
             store.listw = TRUE,
             silent = TRUE)
    {
        if (inherits(xyORdist, "dist")) {
            matdist <- xyORdist
            xy <- dudi.pco(matdist, scannf = FALSE, nf = 2)$li
            if (ncol(xy) == 1)
                xy <- cbind(xy, rep(0, nrow(xy)))
        }
        else {
            matdist <- dist(xyORdist)
            xy <- xyORdist
        }
        
        xy <- as.matrix(xy)
        epsilon <- 10e-6
        
        a <- system.time({
            if (is.null(thresh)) {
                threshh <- give.thresh(matdist)
                if (!silent)
                    cat("Truncation level =", threshh + epsilon, '\n')
            } else {
                threshh <- thresh
                if (!silent)
                    cat("User-provided truncation threshold =",
                        thresh,
                        '\n')
            }
            
            ## compute spatial weighting matrix and associated MEMs
            nb <- dnearneigh(as.matrix(xy), 0, (threshh + epsilon))
            
            spwt <-
                lapply(nbdists(nb, xy), function(x)
                    as.matrix(1 - (x / (4 * threshh)) ^ 2))
            lw <-
                nb2listw(nb,
                         style = "B",
                         glist = spwt,
                         zero.policy = TRUE)
            
            res <-
                scores.listw(lw, MEM.autocor = MEM.autocor, store.listw = store.listw)
        })
        
        a[3] <- sprintf("%2f", a[3])
        if (!silent)
            cat("Time to compute PCNMs =", a[3], " sec", '\n')
        
        attr(res, "call") <- match.call()
        
        return(res)
    }
