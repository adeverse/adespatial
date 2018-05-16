#' Function to create a list of spatial weighting matrices
#' 
#' This function is a user-friendly way to create a list of spatial weighting
#' matrices (SWM) by selecting a set of predefined connectivity and
#' weighting matrices (B and A matrices, respectively). The list can then be fed
#' to the function \code{listw.select} to optimize the selection of the SWM
#' and select the best eigenvector subset within this matrix while controlling
#' the type I error rate.
#' 
#' @details The function allows to construct SWMs based on any combination
#'   of B and A matrices. The B matrices are either graph-based or
#'   distance-based. The function proposes the Delaunay triangulation, Gabriel
#'   graph, relative neighbourhood graph, and the minimum spanning tree criteria
#'   to build a graph-based B matrix. Distance-based SWMs can be built
#'   with the principal coordinates of neighbour matrices (PCNM; Borcard and
#'   Legendre 2002) criteria (see details below), or using another threshold 
#'   distance to define the connected site pairs. The A matrix can be based on a
#'   binary, linear, concave-down, or concave-up function. The linear,
#'   concave-down, and concave-up weighting functions are defined by \eqn{1 -
#'   (D/dmax)}, \eqn{1 - (D/dmax)^y}, and \eqn{1 / D^y}, respectively, where
#'   \code{D} is the euclidean distance between the two sites considered,
#'   \code{dmax} is the maximum euclidean distance between two sites, and 
#'   \code{y} is a user-defined parametre that can either be a single value or a
#'   vector of values. The choice \code{nb = "pcnm"} consists in constructing a
#'   distance-based SWM based on the largest edge of the minimum spanning
#'   tree as a connectivity distance threshold, and then by weighting the links
#'   by the function \eqn{1-(D/(4*t))^2}, where \code{D} is the euclidean
#'   distance between the sites, and \code{t} is the distance threshold below
#'   which two sites are considered connected (Dray et al. 2006). As optimizing
#'   the choice of a SWM has to be done with a p-value correction depending
#'   on the number of candidate SWMs tested (see function
#'   \code{listw.select}), Bauman et al. (2018) strongly encouraged plotting the
#'   concave-down and concave-up weighting functions with several parametre
#'   values in order to only choose the realistic ones to build the candidate W
#'   matrices (e.g., ranging between 0.1 and 1 for the concave-up function, as
#'   values over 1 would make no ecological sense). First visualizing the 
#'   connectivity schemes with the \code{createlistw} function may also help
#'   choosing the B matrices to select for the \code{listw.candidates} function.
#'   
#' @param coord Vector, matrix, or dataframe of point coordinates
#' @param style Coding scheme style (see \code{nb2listw} of the \code{spdep}
#'   package). Can take values 'W', 'B', 'C', 'U', 'minmax', and 'S'; default is
#'   'B'
#' @param nb Defines how the B matrix (connectivity) is build:
#'  * \code{del} Delaunay triangulation
#'  * \code{gab} Gabriel's graph 
#'  * \code{rel} Relative neighbourhood graph
#'  * \code{mst} Minimum spanning tree
#'  * \code{pcnm} Distance-based SWM based on the principal 
#'   coordinates of neighbour matrices (PCNM) criteria (see
#'   'Details')
#'  * \code{dnear} Distance-based
#'  
#' @param d2 Only considered if \code{nb = "dnear"}. It defines the connectivity
#'   distance threshold below which two sites are connected (i.e., maximum distance between two
#'   neighbors. It can either be a single value or a vector of values, in which case a 
#'   different SWM will be generated for each threshold value. The default value is the 
#'   minimum distance keeping all points connected (i.e., the largest edge of the minimum 
#'   spanning tree)
#' @param d1 Only considered if \code{nb = "dnear"}. A single value defining the distance beyond which 
#'   two sites are connected (i.e., minimum distance between two neighbor sites). The default 
#'   value is 0 (no constraint on the min distance). \code{d1} must be smaller than \code{d2}
#' @param weights Defines how the A matrix (weighths) is build:
#'  * \code{binary} without weights
#'  * \code{flin} Linear weighting function 
#'  * \code{fdown} Concave-down weighting function(see Details below)
#'  * \code{fup} Concave-up weighting function (see Details below)
#'  
#' @md 
#' @param y_fdown Single value or vector of values of the \code{y} parameter
#'   in the concave-down weighting function; default is 5
#' @param y_fup Single value or vector of values of the \code{y} parameter
#'   in the concave-up weighting function; default is 0.5
#'   
#' @return A list of SWMs. Each element of the list was built by
#'   \code{nb2listw} (package \code{spdep}) and therefore is of class
#'   \code{listw} and \code{nb}. The name of each element of the list (SWM)
#'   is composed of the corresponding B and A matrices, followed (if any) by the
#'   \code{y} parameter value of the weighting function.
#'   
#' @author David Bauman (\email{dbauman@@ulb.ac.be} or \email{davbauman@@gmail.com}) and Stéphane Dray
#'   
#' @seealso \code{\link{createlistw}}, \code{\link{listw.select}}
#'   
#' @references Bauman D., Fortin M-J., Drouet T. and Dray S. (2018) Optimizing the choice of 
#' a spatial weighting matrix in eigenvector-based methods. Ecology
#' 
#' Borcard D. and Legendre P. (2002) All-scale spatial analysis of 
#' ecological data by means of principal coordinates of neighbour matrices.
#' Ecological Modelling, 153, 51--68
#'   
#' Dray S., Legendre P. and Peres-Neto P. R. (2006) Spatial modeling: a
#' comprehensive framework for principal coordinate analysis of neighbor
#' matrices (PCNM). Ecological Modelling, 196, 483--493
#'   
#' @keywords spatial
#'   
#' @examples
#' ### Create 100 random sampling locations in a squared grid of 120 x 120:
#' xy <- matrix(nrow = 100, ncol = 2)
#' xy[, 1] <- sample(c(1:120), 100, replace = FALSE)
#' xy[, 2] <- sample(c(1:120), 100, replace = FALSE)
#' ### The function listw.candidates is used to build the spatial weighting matrices that
#' ### we want to test and compare (with the listw.select function). We test a Gabriel's graph, 
#' ### a minimum spanning tree, and a distance-based connectivity defined by a threshold
#' ### distance corresponding to the smallest distance keeping all sites connected (i.e., 
#' ### the defaut value of d2). These connectivity matrices are then either not weighted 
#' ### (binary weighting), or weighted by the linearly decreasing function:
#' candidates <- listw.candidates(coord = xy, nb = c("gab", "mst", "dnear"), 
#'                                weights = c("binary", "flin"))
#' names(candidates)                              
#' plot(candidates[[1]], xy)
#' plot(candidates[[3]], xy)
#' ### Construction of a different list of spatial weighting matrices. This time, the
#' ### connexions are defined by a distance-based criterion based on the same threshold
#' ### value, but the connections are weighted by the concave-down function with a y parameter
#' ### varying between 2 and 5, and a concave-up function with a y parametre of 0.2.
#' candidates2 <- listw.candidates(coord = xy, nb = "dnear", weights = c("fdown", "fup"),
#'                                 y_fdown = 1:5, y_fup = 0.2)
#' ### Number of spatial weighting matrices generated:
#' length(candidates2) 
#' ### A single SWM can also easily be generated with listw.candidates:
#' lw <- listw.candidates(xy, nb = "gab", weights = "bin")
#' plot(lw[[1]], xy)
#' 
#' @importFrom spdep tri2nb nb2listw nbdists graph2nb gabrielneigh relativeneigh dnearneigh
#' @export

"listw.candidates" <- function (coord, 
    style = "B", 
    nb = c("del", "gab", "rel", "mst", "pcnm", "dnear"),
    d1 = 0, 
    d2,
    weights = c("binary", "flin", "fup", "fdown"),
    y_fdown = 5, 
    y_fup = 0.5) {
    
    nb <- match.arg(nb, several.ok = TRUE)
    weights <- match.arg(weights, several.ok = TRUE)
    
    
    if(length(nb) == 0) 
        stop("No connectivity matrix selected")
    if (length(weights) == 0) 
        stop("No weighting matrix selected")
    
    if (anyNA(coord)) stop("NA entries in coord")
    
    coord.mat <- as.matrix(coord)
    xy.d1 <- dist(coord.mat)  
    
    res <- list()
    ## Manage the case of pcnm separately (as it is not concerned by weights)
    if ("pcnm" %in% nb) {
        f <- function (D, t) { 1-(D/(4*t))^2 }           # PCNM criterion
        lowlim <- give.thresh(xy.d1)
        matB <- dnearneigh(x = coord.mat, d1 = 0, d2 = lowlim)
        res[[1]] <- nb2listw(matB, style = style, 
                    glist = lapply(nbdists(matB, coord.mat),
                        f, t = lowlim))
        names(res) <- "DBEM_PCNM"
        nb <- nb[-match("pcnm", nb)]
        if(length(nb) == 0)
            return(res)
    }
    
    
    .addweights <- function(nb.object, nb.dist, coord.mat, weights, style, Prefix,
        y_fdown, y_fup) {
        ## this utility function returns a list of listw objects created 
        ## using a single nb.object and several weights
        res <- list()
        res.names <- c()
        f1 <- function (D, dmax)     { 1 - (D/dmax) }       # Linear function
        f2 <- function (D, dmax, y)  { 1 - (D/dmax)^y }     # Concave-down function
        f3 <- function (D, y)        { 1 / D^y }            # Concave-up function
        
        if ("binary" %in% weights) {
            res <- c(res, list(nb2listw(nb.object, style = style)))
            res.names <- c(res.names, paste(Prefix, "Binary", sep = "_"))
        }
        
        if ("flin" %in% weights) {
            max.del <- max(unlist(nbdists(nb.object, coord.mat)))
            
            res <- c(res, list(nb2listw(nb.object, style = style, 
                glist = lapply(nb.dist, f1, dmax = max.del))))
            res.names <- c(res.names, paste(Prefix, "Linear", sep = "_"))
        }
        
        if ("fdown" %in% weights) {
            max.del <- max(unlist(nbdists(nb.object, coord.mat)))
            for (i in y_fdown) {
                res <- c(res, list(nb2listw(nb.object, style = style, 
                    glist = lapply(nbdists(nb.object, coord.mat), 
                        f2, y = i, dmax = max.del))))
                res.names <- c(res.names, paste(Prefix, "Down", i, sep = "_"))
                
            }
        }
        
        
        if ("fup" %in% weights) {
            for (i in y_fup) {
                res <- c(res, list(nb2listw(nb.object, style = style, 
                    glist = lapply(nbdists(nb.object, coord.mat), 
                        f3, y = i))))
                res.names <- c(res.names, paste(Prefix, "Up", i, sep = "_"))
                
            }
        }
        
        names(res) <- res.names
        return(res)
    }
    
    
    
    if ("del" %in% nb) {
        nb.object <- tri2nb(coord.mat)
        nb.dist <- nbdists(nb.object, coord.mat)
        res <- c(res, .addweights(nb.object, nb.dist, coord.mat, weights, style, 
            Prefix = "Delaunay", y_fdown, y_fup))
    }
    
    
    if ("gab" %in% nb) {
        nb.object <- graph2nb(gabrielneigh(coord.mat, nnmult = 5), sym = TRUE)
        nb.dist <- nbdists(nb.object, coord.mat)
        res <- c(res, .addweights(nb.object, nb.dist, coord.mat, weights, style, 
            Prefix = "Gabriel", y_fdown, y_fup))
    }
    
    
    if ("rel" %in% nb) {
        nb.object <- graph2nb(relativeneigh(coord.mat, nnmult = 5), sym = TRUE)
        nb.dist <- nbdists(nb.object, coord.mat)
        res <- c(res, .addweights(nb.object, nb.dist, coord.mat, weights, style, 
            Prefix = "Relative", y_fdown, y_fup))
    }
    
    if ("mst" %in% nb) {
        nb.object <- mst.nb(xy.d1)
        nb.dist <- nbdists(nb.object, coord.mat)
        res <- c(res, .addweights(nb.object, nb.dist, coord.mat, weights, style, 
            Prefix = "MST", y_fdown, y_fup))
    }
    
    if ("dnear" %in% nb) {
        # If "dnear" and no value specified for d2, then d2 is set to be the
        # largest edge of the minimum spanning tree (minimum distance keeping all points connected)
        if (missing(d2))
            d2 <- give.thresh(xy.d1)
    
        nb.dnear <- lapply(d2, dnearneigh, x = coord.mat, d1 = d1)
        for(nb.i in 1:length(nb.dnear)){
            nb.object <- nb.dnear[[nb.i]]
            nb.dist <- nbdists(nb.object, coord.mat)
            res <- c(res, .addweights(nb.object, nb.dist, coord.mat, weights, style, 
                Prefix = paste("Dnear", round(d2[nb.i], 2), sep = ""), y_fdown, y_fup))
            
            
        }
        
    }
    
    return(res)
    
}
