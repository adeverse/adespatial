##
## /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
## |                                                            |
## |  CONSTRAINED HIERARCHICAL CLUSTERING                       |
## |  using (user-specified) criterion                          |
## |                                                            |
## |  C implementation of the Lance and Williams (1967) method  |
## |  for hierarchical clustering with or without the spatial   |
## |  contiguity constraint.                                    |
## |                                                            |
## |  Guillaume Guénard, Université de Montréal, Québec, Canada |
## |  August 2018 - February 2020                               |
## |                                                            |
## \-----------------------------------------------------------*/
##
#' Space- And Time-Constrained Clustering
#'
#' Function \code{constr.hclust} carries out space-constrained or
#' time-constrained agglomerative clustering from a multivariate dissimilarity
#' matrix.
#'
#' @param d A \code{\link{dist}-class} dissimilarity (distance) matrix
#' @param method The agglomeration method to be used (default: "ward.D2"; see
#' details)
#' @param links A list of edges (or links) connecting the points. May be omitted
#' in some cases; see details and examples
#' @param coords Coordinates of the observations (data rows) in \code{d} (for
#' data plotting purposes; may be omitted: see details and examples)
#' @param beta The beta parameter for beta-flexible clustering (default:
#' \code{beta = -0.25})
#' @param chron Logical (TRUE or FALSE) indicating whether a chronological (i.e.
#' time-constrained) clustering should be calculated (default:
#' \code{chron = FALSE})
#'
#' @return A \code{\link{constr.hclust-class}} object.
#' 
#' @details The agglomeration method to be used should be (an unambiguous
#' abbreviation of) one of \code{"ward.D"},
#' \code{"ward.D2"}, \code{"single"}, \code{"complete"},
#' \code{"average"} (UPGMA), \code{"mcquitty"} (WPGMA),
#' \code{"centroid"} (UPGMC), \code{"median"} (WPGMC), or
#' \code{"flexible"}. Method \code{"ward.D2"} (default) implements the
#' Ward (1963) clustering criterion, method
#' \code{"ward.D"} does not (Murtagh and Legendre, 2014).
#'
#' Agglomerative clustering can be carried out with a constraint of spatial or
#' temporal contiguity. This means that only the objects that are linked in
#' \code{links} are considered to be candidates for clustering: the next pair of
#' objects to cluster will be the pair that has the lowest dissimilarity value
#' among the pairs that are linked.
#'
#' The same rule applies during the subsequent clustering steps, which involve
#' groups of objects: the list of links is updated after each agglomeration
#' step. All objects that are neighbours of one of the components that have
#' fused are now neighbours of the newly formed cluster.
#'
#' The edges (links) are specified using argument \code{links}, which can be an
#' object of class \code{nb} (see, e.g., \code{\link{tri2nb}}), an object of
#' class \code{listw} (see, e.g., \code{\link{nb2listw}}), a two-element
#' \code{list} or an object coercible as a such (e.g., a two-column
#' \code{dataframe}), or a two-column matrix with each row representing an edge
#' and the columns representing the two ends of the edges. For lists with more
#' than two elements, as well as dataframes or matrices with more than
#' two-columns only the first two elements or columns are used for the analysis.
#' The edges are interpreted as being non directional; there is no need to
#' specify an edge going from point a to point b and one going from point b to
#' point a. While doing so is generally inconsequential for the analysis, it
#' carries some penality in terms of computation time. It is a good practice to
#' place the nodes in increasing order of numbers from the top to the bottom and
#' from the left to the right of the list but this is not mandatory. A word of
#' caution: in cases where clusters with identical minimum distances occur, the
#' order of the edges in the list may have an influence on the result.
#' Alternative results would be statistically equivalent.
#' 
#' When argument \code{link} is omitted, regular (unconstrained) clustering is
#' performed and a \code{\link{hclust}-class} object is returned unless
#' argument \code{chron = TRUE}. When argument \code{chron = TRUE},
#' chronological clustring is performed, taking the order of observations as
#' their positions in the sequence. Argument \code{links} is not used when
#' \code{chron = TRUE}. Argument \code{chron} allows one to perform a
#' chronological clustring in the case where observations are ordered
#' chronologically. Here, the term "chronologically" should not be taken
#' restrictively: the method remains applicable to other sequential data sets
#' such as spatial series made of observations along a transect.
#' 
#' Memory storage and time to compute constrained clustering for N objects. ---
#' The Lance and Williams algorithm for agglomerative clustering uses
#' dissimilarity matrices. The amount of memory needed to store the distances
#' among N observations as 64-bit double precision floating point variables
#' (IEEE 754) is 8*N*(N-1)/2 bytes. For example, a dissimilarity matrix among
#' 22500 observations would require 2 024 910 000 bytes (1.89 GiB) of storage
#' whereas one among 100 000 observations would take up 39 999 600 000 bytes
#' (37.25 GiB). The implementation in this function needs to cache a copy of the
#' dissimilarity matrix as its elements are modified following each merging of
#' the closest clusters or singletons, thereby doubling the amounts of required memory
#' shown above. Memory needed to store the other informations associated with
#' the clustering is much smaller. Users should make sure to have
#' the necessary memory space (and system stability) before attempting to
#' analyze large data sets. What is considered a large amount of memory has
#' increased over time as computer hardware evolved with time. We let users
#' apply contemporary common sense as to what sample sizes represent manageable
#' clustering problems. Computation time grows with N at roughly the same speed
#' as memory storage requirement to store the dissimilarity matrices increases. 
#' See the Benchmarking example below.
#' 
#' With large data sets, a manageable output describing the classification of
#' the sites is obtained with function \code{\link{cutree}}(x, k) where k is the
#' number of groups.
#' 
#' @author Pierre Legendre \email{pierre.legendre@@umontreal.ca} 
#' (preliminary version coded in R)
#' and Guillaume Guénard \email{guillaume.guenard@@umontreal.ca} 
#' (present version mostly coded in C)
#'
#' @seealso \code{\link{plot.constr.hclust}}, \code{\link{hclust}},
#' \code{\link{cutree}}, and \code{\link{ScotchWhiskey}}
#' 
#' @references
#' Legendre, P. and L. Legendre. 2012. Numerical ecology, 3rd English edition.
#' Elsevier Science BV, Amsterdam.
#' 
#' Murtagh, F. and P. Legendre. 2014. Ward’s hierarchical agglomerative
#' clustering method: which algorithms implement Ward’s criterion? Journal of
#' Classification 31: 274-295. doi: 10.1007/s00357-014-9161-z
#' 
#' Ward, J. H. 1963. Hierarchical grouping to optimize an objective function.
#' Journal of the American Statistical Association 58: 236-244.
#'
#' @examples
#' ##
#' ### First example: Artificial map data from Legendre & Legendre
#' ###                (2012, Fig. 13.26): n = 16
#' ##
#' dat <- c(41,42,25,38,50,30,41,43,43,41,30,50,38,25,42,41)
#' coord.dat <- matrix(c(1,3,5,7,2,4,6,8,1,3,5,7,2,4,6,8,
#'                       4.4,4.4,4.4,4.4,3.3,3.3,3.3,3.3,
#'                       2.2,2.2,2.2,2.2,1.1,1.1,1.1,1.1),16,2)
#' ##
#' ### Obtaining a list of neighbours:
#' library(spdep)
#' listW <- nb2listw(tri2nb(coord.dat), style="B")
#' links.mat.dat <- listw2mat(listW)
#' neighbors <- listw2sn(listW)[,1:2]
#' ##
#' ### Calculating the (Euclidean) distance between points:
#' D.dat <- dist(dat)
#' ##
#' ### Display the points:
#' plot(coord.dat, type='n',asp=1)
#' title("Delaunay triangulation")
#' text(coord.dat, labels=as.character(as.matrix(dat)), pos=3)
#' for(i in 1:nrow(neighbors))
#'     lines(rbind(coord.dat[neighbors[i,1],],
#'           coord.dat[neighbors[i,2],]))
#' ##
#' ### Unconstrained clustring by hclust:
#' grpWD2_hclust <- hclust(D.dat, method="ward.D2")
#' plot(grpWD2_hclust, hang=-1)
#' ##
#' ### Clustering without a contiguity constraint;
#' ### the result is represented as a dendrogram:
#' grpWD2_constr_hclust <- constr.hclust(D.dat, method="ward.D2")
#' plot(grpWD2_constr_hclust, hang=-1)
#' ##
#' ### Clustering with a contiguity constraint described by a list of
#' ### links:
#' grpWD2cst_constr_hclust <-
#'     constr.hclust(
#'         D.dat, method="ward.D2",
#'         neighbors, coord.dat)
#' ##
#' ### To visualize using hclust's plotting method:
#' ### stats:::plot.hclust(grpWD2cst_constr_hclust, hang=-1)
#' ##
#' ### Plot the results on a map with k=3 clusters:
#' plot(grpWD2cst_constr_hclust, k=3, links=TRUE, las=1, xlab="Eastings",
#'      ylab="Northings", pch=21L, cex=3, lwd=3)
#' ##
#' ### Generic functions from hclust can be used, for instance to obtain
#' ### a list of members of each cluster:
#' cutree(grpWD2cst_constr_hclust, k=3)
#' ##
#' ### Now with k=5 clusters:
#' plot(grpWD2cst_constr_hclust, k=5, links=TRUE, las=1, xlab="Eastings",
#'      ylab="Northings", pch=21L, cex=3, lwd=3)
#' cutree(grpWD2cst_constr_hclust, k=5)
#' ##
#' ## End of the artificial map example
#' 
#'
#' ### Second example: Fish community composition along the Doubs River,
#' ### France. The sequence is analyzed as a case of chronological
#' ### clustering, substituting space for time.
#' ##
#' data(doubs, package="ade4")
#' Doubs.D <- dist.ldc(doubs$fish, method="hellinger")
#' grpWD2cst_fish <- constr.hclust(Doubs.D, method="ward.D2", chron=TRUE,
#'                                 coords=as.matrix(doubs$xy))
#' plot(grpWD2cst_fish, k=5, las=1, xlab="Eastings (km)",
#'      ylab="Northings (km)", pch=21L, cex=3, lwd=3)
#' ##
#' ### Repeat the plot with other values of k (number of groups)
#' ##
#' ## End of the Doubs River fish assemblages example
#' 
#'
#' ### Third example: Scotch Whiskey distilleries clustered using tasting scores
#' ### (nose, body, palate, finish, and the four distances combined) constrained
#' ### with respect to the distillery locations in Scotland.
#' ##
#' ## Documentation file about the Scotch Whiskey data: ?ScotchWhiskey
#' ##
#' data(ScotchWhiskey)
#' ### Cluster analyses for the nose, body, palate, and finish D
#' ### matrices:
#' grpWD2cst_ScotchWhiskey <-
#'     lapply(
#'         ScotchWhiskey$dist,    ## A list of distance matrices
#'         constr.hclust,         ## The function called by function lapply
#'         links=ScotchWhiskey$neighbors@data,         ## The list of links
#'         coords=ScotchWhiskey$geo@coords/1000
#'     )
#' ##
#' ### The four D matrices (nose, body, palate, finish), represented as
#' ### vectors in the ScotchWiskey data file, are combined as follows to
#' ### produce a single distance matrix integrating all four types of
#' ### tastes:
#' Dmat <- ScotchWhiskey$dist
#' ScotchWhiskey[["norm"]] <-
#'     sqrt(Dmat$nose^2 + Dmat$body^2 + Dmat$palate^2 + Dmat$finish^2)
#' ##
#' ### This example shows how to apply const.clust to a single D matrix when
#' ### the data file contains several matrices.
#' grpWD2cst_ScotchWhiskey[["norm"]] <-
#'     constr.hclust(
#'         d=ScotchWhiskey[["norm"]],method="ward.D2",
#'         ScotchWhiskey$neighbors@data,
#'         coords=ScotchWhiskey$geo@coords/1000
#'     )
#' ##
#' ### A fonction to plot the Whiskey clustering results
#' plotWhiskey <- function(wh, k) {
#'    par(fig=c(0,1,0,1))
#'    plot(grpWD2cst_ScotchWhiskey[[wh]], k=k, links=TRUE, las=1,
#'         xlab="Eastings (km)", ylab="Northings (km)", cex=0.1, lwd=3,
#'         main=sprintf("Feature: %s",wh))
#'    text(ScotchWhiskey$geo@coords/1000,labels=1:length(ScotchWhiskey$geo))
#'    legend(x=375, y=700, lty=1L, lwd=3, col=rainbow(1.2*k)[1L:k],
#'           legend=sprintf("Group %d",1:k), cex=1.25)
#'    SpeyZoom <- list(xlim=c(314.7,342.2), ylim=c(834.3,860.0))
#'    rect(xleft=SpeyZoom$xlim[1L], ybottom=SpeyZoom$ylim[1L],col="#E6E6E680",
#'         xright=SpeyZoom$xlim[2L], ytop=SpeyZoom$ylim[2L], lwd=2, lty=1L)
#'    par(fig=c(0.01,0.50,0.46,0.99), new=TRUE)
#'    plot(grpWD2cst_ScotchWhiskey[[wh]], xlim=SpeyZoom$xlim,
#'         ylim=SpeyZoom$ylim, k=k, links=TRUE, las=1, xlab="", ylab="",
#'         cex=0.1, lwd=3, axes=FALSE)
#'    text(ScotchWhiskey$geo@coords/1000,labels=1:length(ScotchWhiskey$geo))
#'    rect(xleft=SpeyZoom$xlim[1L], ybottom=SpeyZoom$ylim[1L],
#'         xright=SpeyZoom$xlim[2L], ytop=SpeyZoom$ylim[2L], lwd=2, lty=1L)
#' }
#' ##
#' ### Plot the clustering results on the map of Scotland for 5 groups.
#' ### The inset map shows the Speyside distilleries in detail:
#' plotWhiskey("nose", 5L)
#' plotWhiskey("body", 5L)
#' plotWhiskey("palate", 5L)
#' plotWhiskey("finish", 5L)
#' plotWhiskey("norm", 5L)
#' ##
#' ## End of the Scotch Whiskey tasting data example
#'
#'
#' \dontrun{
#' ### Benchmarking example
#' ### Benchmarking can be used to estimate computation time for different
#' ### values of N. 
#' ### Computing time grows with N at roughly the same speed as the memory
#' ### storage requirements to store the dissimilarity matrices.
#' ##
#' require(magrittr)
#' require(pryr)
#' ##
#' benchmark <- function(nobj) {
#'     # Argument -
#'     # nobj : Number of objects in simulation runs
#'     res <- matrix(NA,length(nobj),3) %>% as.data.frame
#'     colnames(res) <- c("N.objects","Storage (MiB)","Time (sec)")
#'     res[,1L] <- nobj
#'     ## resources <- list()
#'     for(i in 1:length(nobj)) { 
#'         N <- nobj[i]
#'         coords.mem <- cbind(x=runif(N,-1,1),y=runif(N,-1,1))
#'         dat.mem <- runif(N,0,1)
#'         if(i>1L) rm(D.mem) ; gc()
#'         D.mem <- try(dat.mem %>% dist)  #; gc()
#'         if(any(class(D.mem)=="try-error"))
#'             break
#'         neighbors.mem <-
#'             (coords.mem %>%
#'                  tri2nb %>%
#'                  nb2listw(style="B") %>%
#'                  listw2sn)[,1:2]
#'         {start.time = Sys.time()
#'             res.mem <- try(constr.hclust(D.mem, method="ward.D2",
#'                                          neighbors.mem))
#'             end.time = Sys.time()}
#'         if(any(class(res.mem)=="try-error"))
#'             break
#'         res[i,2L] <- (2*object_size(D.mem) + object_size(neighbors.mem) +
#'                           object_size(res.mem))/1048576  # n. bytes per MiB
#'         res[i,3L] <- end.time-start.time
#'     }
#'     res[["N.objects"]] <- as.integer(res[["N.objects"]])
#'     res
#' }
#' res <- benchmark(nobj=c(1000,2000,5000,10000,20000,50000,100000))
#' ##
#' ### Plotting the results:
#' ok <- res %>% apply(1L, function(x) !x %>% is.na %>% any)
#' par(mar=c(3,6,2,2),mfrow=c(2L,1L))
#' barplot(height = res[ok,"Time (sec)"], names.arg= res[ok,"N.objects"],
#'         ylab="Time (seconds)\n",xlab="",las=1L,log="y")
#' par(mar=c(5,6,0,2))
#' barplot(height = res[ok,"Storage (MiB)"], names.arg= res[ok,"N.objects"],
#'         ylab="Total storage (MB)\n",xlab="Number of observations",
#'         las=1L,log="y")
#' ##
#' ### Examine the output file
#' res
#' }
#' ## End of the benchmarking example
#' ##
#' ### End of examples
#' 
#' @useDynLib adespatial, .registration = TRUE
#' @importFrom graphics segments points
#' @importFrom grDevices dev.hold dev.flush rainbow
#' @importFrom spdep nb2listw tri2nb listw2mat listw2sn
#' @importFrom stats dist hclust cutree
#' 
#' @export constr.hclust
constr.hclust <- function(d, method = "ward.D2", links, coords, beta = -0.25,
                          chron = FALSE) {
    METHODS <- c("ward.D", "ward.D2", "single", "complete", "average",
                 "mcquitty", "centroid", "median", "flexible")
    i.meth <- pmatch(method, METHODS)
    if (is.na(i.meth))
        stop("invalid clustering method", paste("", method))
    if (i.meth == -1)
        stop("ambiguous clustering method", paste("", method))
    n <- as.integer(attr(d, "Size"))
    if (is.null(n))
        stop("invalid dissimilarities")
    if (n < 2)
        stop("must have n >= 2 objects to cluster")
    len <- as.integer(n*(n - 1)/2)
    if (length(d) != len)
        (if (length(d) < len)
             stop
         else warning)("dissimilarities of improper length")
    if ((beta < -1) | (beta >= 1)) stop("Flexible clustering: (-1 <= beta < 1)")
    storage.mode(d) <- "double"
    if (missing(links)) {
        nl <- 0L
        links <- integer(0L)
        type <- if (chron) 2L else 0L
    } else {
        if(class(links)[1L]=="nb")
            links <- nb2listw(links, style="B")
        if(class(links)[1L]=="listW")
            links <- listw2sn(links)[1L:2L]
        if(is.list(links)||is.data.frame(links)) {
            nl <- length(links[[1L]])
            links <- unlist(links)
            storage.mode(links) <- "integer"
        } else if(is.matrix(links)) {
            if(ncol(links) < 2L)
                stop("When 'links' is a matrix, it must have 2 or more columns")
            nl <- nrow(links)
            links <- as.integer(links)
        } else
            stop("'links' must be a list, dataframe, or matrix")
        type <- 1L
    }
    ## This handling of the parameters is only valid for the flexible
    ## clustering. To implement other methods involving (a) parameter value(s),
    ## replace the next statement with one suitable for every cases where (a)
    ## parameter(s) (is/)are involved.
    pars <- as.double(c((1-beta)/2,beta))
    ##
    hcl <- .C("cclust",
              as.integer(n),             ## n
              merge = integer(2*(n-1L)),
              height = double(n-1L),
              order = integer(n),
              d,                         ## diss0
              as.integer(nl),            ## nl
              links,                     ## linkl
              as.integer(i.meth),        ## method
              pars,                      ## alpha and beta
              as.integer(type)           ## type
              )[2L:4L]
    dim(hcl$merge) <- c((n-1L),2L)
    if(missing(coords)) {
        if(chron) {
            coords <- cbind(x=0:(n-1),y=0)
        } else coords <- NULL
    } else {
        if(NCOL(coords) < 2L) {
            coords <- cbind(x=coords,y=0)
        } else coords <- as.matrix(coords)
    }
    return(
        structure(
            c(hcl,
              list(labels = attr(d, "Labels"),
                   method = METHODS[i.meth], call = match.call(),
                   dist.method = attr(d, "method"),
                   links = if(length(links)) matrix(links,nl,2L) else NULL,
                   coords = coords)),
            class = c("constr.hclust","hclust")))
}
##
