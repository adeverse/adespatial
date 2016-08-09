#'Combine dbMEM matrices corresponding to groups of sites
#'
#'This function reads a file containing the Cartesian coordinates of
#'  sites forming different groups on the map, and constructs a combined staggered matrix
#'  of dbMEM spatial eigenvectors, ready for use in RDA.
#'  The method was first described and used in Declerck et al. (2011) and summarized in
#'  the Borcard et al. (2011) book, section 7.4.3.5. These publications provided
#'  preliminary versions of the present function. The present version is more completely
#'  documented. Furthermore, it uses the \code{\link{dbmem}} function of the
#'  \code{adespatial} package for computation of the eigenfunctions.
#'
#'@param coord Optional file containing the Cartesian coordinates of the sites.
#'@param D.mat Optional distance matrix provided by user, class \code{matrix} or
#'  \code{dist}. If \code{D.mat=NULL}, the geographic distance matrix will be computed
#'  from the coordinates provided in file \code{coord}.
#'@param nsites A vector containing the number of sites per group.
#'
#'@return A matrix with \code{n} rows containing a set of \code{k} staggered matrices of
#'  dbMEM eigenfunctions in its diagonal portion; \code{n} is the total number of sites
#'  in the study and \code{k} is the number of groups. Each small matrix contains
#'  the dbMEM functions, modelling positive spatial correlation, describing the spatial
#'  relationships among the sites of a group. The remainder of the matrix is filled with
#'  zeros. Zero is the mean value of all eigenfunctions describing within-group
#'  relationships. This means that during the calculation of RDA, the sites of a focus
#'  group will have, with each other, relationships described by the dbMEM eigenfunctions
#'  of that group, whereas the sites outside that group will have weights of 0 in the
#'  regressions that concern these eigenfunctions.
#'
#'@details
#'The geographic positions of the sites are provided either in a file of geographic
#'  coordinates \code{coord} or as a geographic distance matrix \code{D.mat}.
#'
#'The sites must, of course, be in the same order in file \code{coord} (or in file
#'  \code{D.mat}) and in the response data file used in the RDA. All sites of a group must
#'  be together in these two files, i.e. not interspersed. The numbers of sites in the
#'  groups are provided in vector \code{nsites}. See example.
#'
#'File vector \code{coord}, if provided, must contain Cartesian coordinates of the sites,
#'  not coordinates in degrees. The Euclidean distance computed from the geographic
#'  coordinates is a meaningful representation of the geographic relationships only if the
#'  coordinates are Cartesian. Geodetic Cartesian coordinates can be derived from Lat-Lon
#'  data in degrees using the function \code{geoXY} of the \code{SoDA} package. Beware of
#'  UTM coordinates if the sites are not all located in the same UTM zone; UTM coordinates
#'  are Cartesian only within an UTM zone. See
#'  https://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system.
#'
#'@author  Pierre Legendre \email{pierre.legendre@@umontreal.ca}, 2010. Adaptation to adespatial: Daniel Borcard and Pierre
#'  Legendre, 2016
#'
#'@references
#'
#'Borcard, D., F. Gillet and P. Legendre. 2011. Numerical ecology with R. Use R! series,
#'  Springer Science, New York.
#'
#'Declerck, S. A. J., J. S. Coronel, P. Legendre & L. Brendonck. 2011. Scale dependency of
#'  processes structuring metacommunities of cladocerans in temporary pools of High-Andes
#'  wetlands. Ecography 34: 296-305.
#'
#'@seealso \code{\link{dbmem}}
#'
#'@keywords spatial
#'
#'@examples {
#'  # Generate random coordinates for 35 sites forming 6 distinct groups on the map
#'  Easting <- runif(35)+c(rep(0,6),rep(1.5,7),rep(3,6), rep(0,5),rep(1.5,5),rep(3,6))
#'  Northing<- runif(35)+c(rep(2.8,6),rep(2.3,7),rep(2.8,6), rep(0,5),rep(0.5,5),rep(0,6))
#'  cartesian <- cbind(Easting,Northing)
#'  rownames(cartesian) <- paste("S",1:nrow(cartesian),sep='')
#'  nsites.per.group <- c(6,7,6,5,5,6)
#'
#'  result <- create.dbMEM.model(coord=cartesian, nsites=nsites.per.group)
#'
#'  # Draw a map to check the coding of the sites into the groups
#   # First, expand vector 'nsites.per.group' into a 'site.codes' vector with 35 values
#'  site.codes <- unlist(apply(cbind(1:6),1,n=nsites.per.group,function(a,n) rep(a,n[a])))
#'
#'  col.vec <- c("green3","gray99","orange2","gold1","brown3","gray70")
#'  plot(cartesian, pch=22, col="black", bg=col.vec[site.codes], cex=2, ylim=c(0,4),asp=1)
#'  text(cartesian,labels=rownames(cartesian), cex=0.5, pos=3)
#'
#'  # Examine the staggered matrix of dbMEM eigenfunctions
#'  # Not run:
#'  result
#'}
#'
#'
#' @export create.dbMEM.model
#' @importFrom stats as.dist

'create.dbMEM.model' <-
    function(coord = NULL,
        D.mat = NULL,
        nsites)
    {
        if (is.null(coord) & is.null(D.mat))
            stop("Geographic information must be provided in 'coord' or in 'D.mat'")
        #
        if (is.null(D.mat))
            D.mat <- dist(coord)
        D.mat <-
            as.matrix(D.mat) # Necessary step before selection of blocks of distances
        #
        if (!is.null(coord)) {
            n <- nrow(coord)
        } else {
            n <- nrow(D.mat)
        }
        if (sum(nsites) != n)
            stop("Vector nsites does not sum to nrow(coord) or nrow(D.mat)")
        if (min(nsites) == 1)
            stop("At least one group contains a single site")
        #
        out <- matrix(0, n, n)
        end <- 0
        end.mem <- 0
        for (k in 1:length(nsites))
        {
            start <- end + 1
            end <- end + nsites[k]
            tmp <- as.dist(D.mat[start:end, start:end])
            res <-
                dbmem(tmp, MEM.autocor = "positive", store.listw = FALSE)
            dbMEM <- as.matrix(res)
            n.mem <- ncol(dbMEM)
            out[start:end, (end.mem + 1):(end.mem + n.mem)] <- dbMEM
            end.mem <- end.mem + n.mem
        }
        out <- out[, 1:end.mem]
        if (is.null(rownames(coord)))
            rownames(out) <-
            rownames(out, do.NULL = FALSE, prefix = "Site.")
        else
            rownames(out) <- rownames(coord)
        colnames(out) <-
            colnames(out, do.NULL = FALSE, prefix = "dbMEM.")
        out
    }