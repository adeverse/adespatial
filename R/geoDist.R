## **************************************************************************
##
##    This R code file has been repackaged from package SoDA 1.0-6 from the
##    CRAN archives (https://cran.r-project.org/src/contrib/Archive/SoDA/,
##    downloaded 2025-03-17 11:04:00) by:
##
##    Guillaume Guénard <guillaume.guenard@umontreal.ca>
##    Department de sciences biologiques,
##    Université de Montréal
##    Montreal, QC, Canada
##
##    Original authors:
##    John M. Chambers <jmc@r-project.org>
##    2014-06-12
##
## **************************************************************************
##
#' 
#' Geodetic Distances from Latitude and Longitude
#' 
#' Given two sets of points on the earth's surface in latitude, longitude form,
#' returns the geodetic distances in meters between corresponding points.
#' 
#' @name geoDist
#' 
#' @param lat1,lon1,lat2,lon2 Latitude and longitude co-ordinates for the two
#' sets of points.
#' @param NAOK Are \code{NA} values allowed in the co-ordinates?  Default:
#' \code{TRUE}. If so, corresponding elements of the distance will also be
#' \code{NA}.
#' 
#' @details
#' Uses a classic Fortran algorithm implementing a method that allows for the
#' non-spherical shape of the earth.  See comments in the Fortran code for the
#' history of the implementation.
#' 
#' @returns numeric vector of distances, optionally including \code{NA} values
#' if those are allowed and present in any of the coordinates.
#' 
#' @author
#' John Chambers, Stanford Data Science, Stanford University (formerly at Bell
#' Labs). Transfer of this function from SoDA to adespatial by Guillaume
#' Guénard, with permission of John Chambers.
#' 
#' @references
#' Vincenty,T. (1975). Direct and inverse solutions of geodesics on the
#' ellipsoid with application of nested equations. \emph{Survey Review},
#' vol. 23(176):88-94.
#' 
#' @export
"geoDist" <-
function(lat1, lon1, lat2, lon2, NAOK = TRUE) {
  
  n <- unique(c(length(lat1), length(lon1), length(lat2), length(lon2)))
  nok <- n
  
  if(length(n)>1)
    stop("Need all arguments of the same length:  got ",
         paste(n, collapse=", "))
  
  if(n < 1)
    return(numeric())
  
  nas <- is.na(lat1) | is.na(lat2) | is.na(lon1) | is.na(lon2)
  
  if(NAOK) {
    if(any(nas)) {
      ok <- !nas
      lat1 <- lat1[ok]; lon1 <- lon1[ok];
      lat2 <- lat2[ok]; lon2 <- lon2[ok];
      nok <- sum(ok)
    }
  } else if(any(nas))
    stop("NA values found but not allowed")
  
  .Fortran(
    "GEODISTV",
    as.double(lat1),
    as.double(lon1),
    as.double(lat2),
    as.double(lon2),
    dist = double(nok),
    as.integer(nok),
    PACKAGE = "adespatial"
  )$dist -> res
  
  if(NAOK && any(nas)) {
    value <- rep(NA, n)
    value[ok] <- res
    value
  } else
    res
}
