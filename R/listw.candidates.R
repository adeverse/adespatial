#' Function to create a list of spatial weighting matrices
#' 
#' This function is a user-friendly way to create a list of spatial weighting
#' matrices (W matrices) by selecting a set of predefined connectivity and
#' weighting matrices (B and A matrices, respectively). The list can then be fed
#' to the function \code{MEM.modsel} to optimise the selection of the W matrix
#' and select the best eigenvector subset within this matrix while controlling
#' the type I error rate.
#' 
#' @details The function allows to construct W matrices based on any combination
#'   of B and A matrices. The B matrices are either graph-based or
#'   distance-based. The function proposes the Delaunay triangulation, Gabriel
#'   graph, relative neighbourhood graph, and the minimum spanning tree criteria
#'   to build a graph-based B matrix. Distance-based W matrices can be built
#'   with the principal coordinates of neighbour matrices (PCNM; Borcard and
#'   Legendre 2002) criteria (see details below), or using another threshold 
#'   distance to define the connected site pairs. The A matrix can be based on a
#'   binary, linear, concave-down, or concave-up function. The linear,
#'   concave-down, and concave-up weighting functions are defined by \eqn{1 -
#'   (D/dmax)}, \eqn{1 - (D/dmax)^y}, and \eqn{1 / D^y}, respectively, where
#'   \code{D} is the euclidean distance between the two sites considered,
#'   \code{dmax} is the maximum euclidean distance between two sites, and 
#'   \code{y} is a user-defined parametre that can either be a single value or a
#'   vector of values. The argument \code{PCNM} consists in constructing a
#'   distance-based W matrix based on the largest edge of the minimum spanning
#'   tree as a connectivity distance threshold, and then by weighting the links
#'   by the function \eqn{1-(D/(4*t))^2}, where \code{D} is the euclidean
#'   distance between the sites, and \code{t} is the distance threshold below
#'   which two sites are considered connected (Dray et al. 2006). As optimising
#'   the choice of a W matrix has to be done with a p-value correction depending
#'   on the number of W matrix candidates tested (see function
#'   \code{MEM.modsel}), Bauman et al. (2018) strongly encouraged plotting the
#'   concave-down and concave-up weighting functions with several parametre
#'   values in order to only choose the realistic ones to build the W matrix
#'   candidates (e.g., ranging between 0.1 and 1 for the concave-up function, as
#'   values over 1 would make no ecological sense). First visualising the 
#'   connectivity schemes with the \code{createlistw} function may also help
#'   choosing the B matrices to select for the \code{listw.candidates} function.
#'   
#' @param coord Vector, matrix, or dataframe of point coordinates
#' @param style Coding scheme style (see \code{nb2listw} of the \code{spdep}
#'   package). Can take values 'W', 'B', 'C', 'U', 'minmax', and 'S'; default is
#'   'B'
#' @param del Defines whether a B matrix based on the Delaunay triangulation
#'   should be used; default is FALSE. No edge effect correction implemented here (see
#'   discussion of Bauman et al. 2018)
#' @param gab Defines whether a B matrix based on a Gabriel's graph should be
#'   used; default is FALSE
#' @param rel Defines whether a B matrix based on a relative neighbourhood graph
#'   should be used; default is FALSE
#' @param mst Defines whether a B matrix based on a minimum spanning tree should
#'   be used; default is FALSE
#' @param PCNM Defines whether a distance-based W matrix based on the principal 
#'   coordinates of neighbour matrices (PCNM) criteria should be used (see
#'   'Details'); default is FALSE
#' @param DB Defines whether a distance-based W matrix should be built; default
#'   is FALSE
#' @param d2 Only considered if DB is TRUE; defines the connectivity
#'   distance threshold below which two sites are connected (i.e., maximum distance between two
#'   neighbors. It can either be a single value or a vector of values, in which case a 
#'   different W matrix will be generated for each threshold value. The default value is the 
#'   minimum distance keeping all points connected (i.e., the largest edge of the minimum 
#'   spanning tree)
#' @param d1 Only considered if DB is TRUE; single value defining the distance beyond which 
#'   two sites are connected (i.e., minimum distance between two neighbor sites). The default 
#'   value is 0 (no constraint on the min distance). \code{d1} must be smaller than \code{d2}
#' @param binary Defines whether W matrices based on the selected B matrices
#'   should be created without weights on the connexions; default is FALSE
#' @param flin Defines whether the linear weighting function should be used (see
#'   Details below); default is FALSE
#' @param fconcdown Defines whether the concave-down weighting function should
#'   be used (see Details below); default is FALSE
#' @param fconcup Defines whether the concave-up weighting function should be
#'   used (see Details below); default is FALSE
#' @param y_fconcdown Single value or vector of values of the \code{y} parametre
#'   in the concave-down weighting function; default is 5
#' @param y_fconcup Single value or vector of values of the \code{y} parametre
#'   in the concave-up weighting function; default is 0.5
#'   
#' @return A list of W matrices. Each element of the list was built by
#'   \code{nb2listw} (package \code{spdep}) and therefore is of class
#'   \code{listw} and \code{nb}. The name of each element of the list (W matrix)
#'   is composed of the corresponding B and A matrices, followed (if any) by the
#'   \code{y} parametre value of the weighting function.
#'   
#' @author Bauman David \email{dbauman@@ulb.ac.be} or \email{davbauman@@gmail.com}
#'   
#' @seealso \code{\link{createlistw}}, \code{\link{MEM.modsel}}
#'   
#' @references Bauman D., Fortin M-J., Drouet T. and Dray S. (2018) Otimising the choice of 
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
#' ### we want to test and compare (with the MEM.modsel function). We test a Gabriel's graph, 
#' ### a minimum spanning tree, and a distance-based connectivity defined by a threshold
#' ### distance corresponding to the smallest distance keeping all sites connected (i.e., 
#' ### the defaut value of d2). These connectivity matrices are then either not weighted 
#' ### (binary weighting), or weighted by the linearly decreasing function:
#' candidates <- listw.candidates(coord = xy, gab = TRUE, mst = TRUE, DB = TRUE,
#'                                binary = TRUE, flin = TRUE)
#' names(candidates)                              
#' plot(candidates[[1]], xy)
#' plot(candidates[[3]], xy)
#' ### Construction of a different list of spatial weighting matrices. This time, the
#' ### connexions are defined by a distance-based criterion based on the same threshold
#' ### value, but the connections are weighted by the concave-down function with a y parametre
#' ### varying between 2 and 5, and a concave-up function with a y parametre of 0.2.
#' candidates2 <- listw.candidates(coord = xy, DB = TRUE, fconcdown = TRUE, fconcup = TRUE,
#'                                 y_fconcdown = c(1:5), y_fconcup = 0.2)
#' ### Number of spatial weighting matrices generated:
#' length(candidates2) 
#' 
#' @importFrom spdep tri2nb nb2listw nbdists graph2nb gabrielneigh relativeneigh
#'   dnearneigh
#' @export

"listw.candidates" <- function (coord, 
                                style = "B", 
                                del = FALSE, 
                                gab = FALSE, 
                                rel = FALSE, 
                                mst = FALSE, 
                                PCNM = FALSE,  
                                DB = FALSE, 
                                d1 = 0, 
                                d2,
                                binary = FALSE, 
                                flin = FALSE, 
                                fconcdown = FALSE, 
                                fconcup = FALSE, 
                                y_fconcdown = 5, 
                                y_fconcup = 0.5) {
  
  if (PCNM == FALSE) {
    if (any(del, gab, rel, mst, DB) == FALSE) 
      stop("No connectivity matrix selected")
    if (any(binary, flin, fconcdown, fconcup) == FALSE) 
      stop("No weighting matrix selected")
  }
  if (any(is.na(coord))) stop("NA entries in coord")
  # If DB = TRUE and no value was specified for d2, then d2 is set to be the
  # largest edge of the minimum spanning tree (minimum distance keeping all points connected):
  if (DB == TRUE) {
    class <- class(try(is.vector(d2), TRUE))
    if (class == "try-error") d2 <- give.thresh(dist(coord))
  }
  
  if (length(which(c(flin, fconcdown, fconcup) == TRUE)) != 0) weightfun = TRUE
  else weightfun = FALSE
  
  # Definition of the weighting functions:
  if (weightfun == TRUE) {
    f1 <- function (D, dmax)     { 1 - (D/dmax) }       # Linear function
    f2 <- function (D, dmax, y)  { 1 - (D/dmax)^y }     # Concave-down function
    f3 <- function (D, y)        { 1 / D^y }            # Concave-up function
  }
  
  xy.d1 <- dist(coord)
  
  # Total nb of W matrices to be built:
  # ***********************************
  nbB <- length(which(c(del, gab, rel, mst) == TRUE))
  if (DB == TRUE) nbB <- nbB + length(d2)
  nbw <- nbB
  control_BinLin <- FALSE
  control_f2 <- FALSE
  if (length(which(c(binary, flin) == TRUE)) != 0) {
    control_BinLin <- TRUE
    nbw <- nbB * length(which(c(binary, flin) == TRUE))
  }
  if (weightfun == TRUE) {
    if (fconcdown == TRUE) {
      control_f2 <- TRUE
      yf2 <- length(y_fconcdown)
      if (control_BinLin == TRUE) nbw <- nbw + (nbB * yf2) 
      else nbw <- nbB * yf2
    }
    if (fconcup == TRUE) {
      yf3 <- length(y_fconcup)
      if (control_BinLin == TRUE) nbw <- nbw + (nbB * yf3) 
      else if (control_f2 == TRUE) nbw <- nbw + nbB * yf3 else nbw <- nbB * yf3
    }
  }
  if (PCNM == TRUE) nbw <- nbw + 1
  
  # List for the W matrix candidates:
  listwcand <- vector("list", nbw)
  count <- 0
  
  # Construction of the list of W matrix candidates: 
  # ************************************************
  if (del == TRUE) {
    Y.del <- tri2nb(jitter(as.matrix(coord), factor = 0.001))
    if (binary == TRUE) {
      count <- count + 1
      listwcand[[count]] <- nb2listw(Y.del, style = style)
      names(listwcand)[count] <- "Delaunay_Binary"
    }
    if (weightfun == TRUE) {
      max.del <- max(unlist(nbdists(Y.del, as.matrix(coord))))
      if (flin == TRUE) {
        count <- count + 1
        listwcand[[count]] <- nb2listw(Y.del, style = style, 
                                       glist = lapply(nbdists(Y.del, as.matrix(coord)),
                                                      f1, dmax = max.del))
        names(listwcand)[count] <- "Delaunay_Linear"
      } 
      if (fconcdown == TRUE) {      
        count <- count + 1
        for (i in y_fconcdown) {
          listwcand[[count]] <- nb2listw(Y.del, style = style, 
                                         glist = lapply(nbdists(Y.del, as.matrix(coord)), 
                                                        f2, y = i, dmax = max.del))
          names(listwcand)[count] <- paste("Delaunay_Concave down (y = ", i, ")", sep = "")
          if (i != y_fconcdown[length(y_fconcdown)]) count <- count + 1
        }
      }
      if (fconcup == TRUE) {      
        count <- count + 1
        for (i in y_fconcup) {
          listwcand[[count]] <- nb2listw(Y.del, style = style, 
                                         glist = lapply(nbdists(Y.del, as.matrix(coord)), 
                                                        f3, y = i))
          names(listwcand)[count] <- paste("Delaunay_Concave up (y = ", i, ")", sep = "")
          if (i != y_fconcup[length(y_fconcup)]) count <- count + 1
        }
      }
    }
  }
  if (gab == TRUE) {
    Y.gab <- graph2nb(gabrielneigh(as.matrix(coord), nnmult = 5), sym = TRUE)
    if (binary == TRUE) {
      count <- count + 1
      listwcand[[count]] <- nb2listw(Y.gab, style = style)
      names(listwcand)[count] <- "Gabriel_Binary"
    }
    if (weightfun == TRUE) {
      max.gab <- max(unlist(nbdists(Y.gab, as.matrix(coord))))
      if (flin == TRUE) {
        count <- count + 1
        listwcand[[count]] <- nb2listw(Y.gab, style = style, 
                                       glist = lapply(nbdists(Y.gab, as.matrix(coord)),
                                                      f1, dmax = max.gab))
        names(listwcand)[count] <- "Gabriel_Linear"
      } 
      if (fconcdown == TRUE) {      
        count <- count + 1
        for (i in y_fconcdown) {
          listwcand[[count]] <- nb2listw(Y.gab, style = style, 
                                         glist = lapply(nbdists(Y.gab, as.matrix(coord)),
                                                                f2, y = i, dmax = max.gab))
          names(listwcand)[count] <- paste("Gabriel_Concave down (y = ", i, ")", sep = "")
          if (i != y_fconcdown[length(y_fconcdown)]) count <- count + 1
        }
      }
      if (fconcup == TRUE) {      
        count <- count + 1
        for (i in y_fconcup) {
          listwcand[[count]] <- nb2listw(Y.gab, style = style, 
                                         glist = lapply(nbdists(Y.gab, as.matrix(coord)), 
                                                        f3, y = i))
          names(listwcand)[count] <- paste("Gabriel_Concave up (y = ", i, ")", sep = "")
          if (i != y_fconcup[length(y_fconcup)]) count <- count + 1
        }
      }
    }
  }
  if (rel == TRUE) {
    Y.rel <- graph2nb(relativeneigh(as.matrix(coord), nnmult = 5), sym = TRUE)
    if (binary == TRUE) {
      count <- count + 1
      listwcand[[count]] <- nb2listw(Y.rel, style = style)
      names(listwcand)[count] <- "Rel. neighbourhood_Binary"
    }
    if (weightfun == TRUE) {
      max.rel <- max(unlist(nbdists(Y.rel, as.matrix(coord))))
      if (flin == TRUE) {
        count <- count + 1
        listwcand[[count]] <- nb2listw(Y.rel, style = style, 
                                       glist = lapply(nbdists(Y.rel, as.matrix(coord)), 
                                                      f1, dmax = max.rel))
        names(listwcand)[count] <- "Rel. neighbourhood_Linear"
      } 
      if (fconcdown == TRUE) {      
        count <- count + 1
        for (i in y_fconcdown) {
          listwcand[[count]] <- nb2listw(Y.rel, style = style, 
                                         glist = lapply(nbdists(Y.rel, as.matrix(coord)),
                                                        f2, y = i, dmax = max.rel))
          names(listwcand)[count] <- paste("Rel. neighbourhood_Concave down (y = ", i, ")", 
                                           sep = "")
          if (i != y_fconcdown[length(y_fconcdown)]) count <- count + 1
        }
      }
      if (fconcup == TRUE) {      
        count <- count + 1
        for (i in y_fconcup) {
          listwcand[[count]] <- nb2listw(Y.rel, style = style, 
                                         glist = lapply(nbdists(Y.rel, as.matrix(coord)),
                                                        f3, y = i))
          names(listwcand)[count] <- paste("Rel. neighbourhood_Concave up (y = ", i, ")", 
                                           sep = "")
          if (i != y_fconcup[length(y_fconcup)]) count <- count + 1
        }
      }
    }
  }
  if (mst == TRUE) {
    Y.mst <- mst.nb(xy.d1)
    if (binary == TRUE) {
      count <- count + 1
      listwcand[[count]] <- nb2listw(Y.mst, style = style)
      names(listwcand)[count] <- "Min. spanning tree_Binary"
    }
    if (weightfun == TRUE) {
      max.mst <- max(unlist(nbdists(Y.mst, as.matrix(coord))))
      if (flin == TRUE) {
        count <- count + 1
        listwcand[[count]] <- nb2listw(Y.mst, style = style, 
                                       glist = lapply(nbdists(Y.mst, as.matrix(coord)),
                                                      f1, dmax = max.mst))
        names(listwcand)[count] <- "Min. spanning tree_Linear"
      } 
      if (fconcdown == TRUE) {      
        count <- count + 1
        for (i in y_fconcdown) {
          listwcand[[count]] <- nb2listw(Y.mst, style = style, 
                                         glist = lapply(nbdists(Y.mst, as.matrix(coord)),
                                                                f2, y = i, dmax = max.mst))
          names(listwcand)[count] <- paste("Min. spanning tree_Concave down (y = ", i, ")", 
                                           sep = "")
          if (i != y_fconcdown[length(y_fconcdown)]) count <- count + 1
        }
      }
      if (fconcup == TRUE) {      
        count <- count + 1
        for (i in y_fconcup) {
          listwcand[[count]] <- nb2listw(Y.mst, style = style, 
                                         glist = lapply(nbdists(Y.mst, as.matrix(coord)),
                                                        f3, y = i))
          names(listwcand)[count] <- paste("Min. spanning tree_Concave up (y = ", i, ")", 
                                           sep = "")
          if (i != y_fconcup[length(y_fconcup)]) count <- count + 1
        }
      }
    }
  }
  if (PCNM == TRUE) {
    count <- count + 1
    f <- function (D, t) { 1-(D/(4*t))^2 }           # PCNM criterion
    lowlim <- give.thresh(xy.d1)
    matB <- dnearneigh(lowlim, x = as.matrix(coord), d1 = d1)
    listwcand[[count]] <- nb2listw(matB, style = style, 
                                   glist = lapply(nbdists(matB, as.matrix(coord)),
                                                  f, t = lowlim))
    names(listwcand)[count] <- "DBMEM_PCNM"
  }
  if (DB == TRUE) {
    Y.listDB <- lapply(d2, dnearneigh, x = as.matrix(coord), d1 = d1)
    if (binary == TRUE) {
      count <- count + 1
      for (i in 1:length(d2)) {
        listwcand[[count]] <- nb2listw(Y.listDB[[i]], style = style)
        names(listwcand)[count] <- paste("DB", i, "_Binary", sep = "")
        if (i != length(d2)) count <- count + 1
      }
    }
    if (weightfun == TRUE) {
      nbdist <- lapply(Y.listDB, coords = as.matrix(coord), nbdists)
      unlist <- lapply(nbdist, unlist)
      max.list <- lapply(unlist, max)
      if (flin == TRUE) {
        count <- count + 1
        for (i in 1:length(d2)) {
          listwcand[[count]] <- nb2listw(Y.listDB[[i]], style = style, 
                                         glist = lapply(nbdists(Y.listDB[[i]], 
                                                                as.matrix(coord)),
                                                        f1, dmax = max.list[[i]]))
          names(listwcand)[count] <- paste("DB", i, "_Linear", sep = "")
          if (i != length(d2)) count <- count + 1
        }
      } 
      if (fconcdown == TRUE) { 
        for (j in 1:length(d2)) {
          count <- count + 1
          for (i in y_fconcdown) {
            listwcand[[count]] <- nb2listw(Y.listDB[[j]], style = style, 
                                           glist = lapply(nbdists(Y.listDB[[j]], 
                                                                  as.matrix(coord)),
                                                          f2, y = i, dmax = max.list[[j]]))
            names(listwcand)[count] <- paste("DB", j, "_Concave down (y = ", i, ")", sep = "")
            if (i != y_fconcdown[length(y_fconcdown)]) count <- count + 1
          }
        }
      }
      if (fconcup == TRUE) { 
        for (j in 1:length(d2)) {
          count <- count + 1
          for (i in y_fconcup) {
            listwcand[[count]] <- nb2listw(Y.listDB[[j]], style = style, 
                                           glist = lapply(nbdists(Y.listDB[[j]], 
                                                                  as.matrix(coord)), f3, 
                                                          y = i))
            names(listwcand)[count] <- paste("DB", j, "_Concave up (y = ", i, ")", sep = "")
            if (i != y_fconcup[length(y_fconcup)]) count <- count + 1
          }
        }
      }
    }
  }
  return(listwcand)
}
