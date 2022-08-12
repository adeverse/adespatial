## **************************************************************************
##
##    (c) 2018-2022 Guillaume Guénard
##        Department de sciences biologiques
##        Université de Montréal
##        Montreal, QC, Canada
##
## **************************************************************************
##
#' Scotch Whiskey Data Set
#'
#' Single Malt Scotch whiskeys from 109 distilleries
#'
#' @docType data
#' 
#' @keywords Scotch Whiskey
#' 
#' @name ScotchWhiskey
#' 
#' @usage data(ScotchWhiskey)
#' 
#' @format A list with 12 members:
#' \describe{
#' \item{geo}{A \code{\link{SpatialPointsDataFrame-class}} object containing the
#' geographic coordinates and other information about the distilleries.}
#' \item{colour}{The whiskey colour coded as a 14-level factor.}
#' \item{nose}{A set of 12 nasal notes (boolean).}
#' \item{body}{A set of 8 body notes (boolean).}
#' \item{palate}{A set of 15 palatine notes (boolean).}
#' \item{finish}{A set of 19 finish (or after-taste) notes (boolean).}
#' \item{nbChar}{Number of characteristics attributed to each distillery for
#' each of the four sets of boolean features: nose, body, palate, finish.}
#' \item{listW}{A \code{listw} object (see \code{\link{nb2listw}})
#' containing information about the spatial edges (neighbour links) between the distilleries.}
#' \item{links.mat}{A binary square matrix of the spatial connexions between the
#' distilleries (contiguity matrix).}
#' \item{neighbors}{A \code{\link{SpatialLinesDataFrame-class}} object
#' containing geographic information about the spatial links between the
#' distilleries.}
#' \item{dist}{A list of distance matrices obtained for each of the four sets of
#' boolean features.} }
#' 
#' @details There are 5 data sets: color, nose, body, palate, and finish. The
#' binary (0,1) descriptors are in the same order as on p. 239 of the whisky paper.
#'
#' There are two whiskies in the classification from the Springbank distillery.
#' One pertains to the Islay group, the other to the Western group.
#' 
#' Please let us know of the analyses you have performed with the whiskey data,
#' especially if you intend to publish them.
#'
#' The distance matrices were calculated separately as follows for each tasting
#' data set:
#'
#' D = (1 - S4)^0.5,
#'
#' where S4 is the Simple matching coefficient of Sokal & Michener (1958). This
#' coefficient was called S4 in the Gower & Legendre (1986) paper and S1 in the
#' Legendre & Legendre (2012) book. In package ade4, coefficient
#' D = sqrt(1 - S4) is computed by function \code{\link{dist.binary}} using
#' argument \code{"method=2"}.
#' 
#' @source Pierre Legendre <pierre.legendre@@umontreal.ca> and François-Joseph
#' Lapoints <francois-joseph.lapointe@@umontreal.ca>, Département de sciences
#' biologiques, Université de Montréal, Montréal, Québec, Canada.
#' 
#' @references Lapointe, F.-J. and P. Legendre. 1994. A classification of pure
#' malt Scotch whiskies. Applied Statistics 43: 237-257
#' <http://www.dcs.ed.ac.uk/home/jhb/whisky/lapointe/text.html>.
#'
#' Gower, J.C. and Legendre, P. 1986. Metric and Euclidean properties of
#' dissimilarity coefficients. Journal of Classification, 3, 5-48.
#' 
#' Legendre, P. and Legendre, L. 2012. Numerical Ecology. 3rd English edition.
#' Elsevier Science BV, Amsterdam.
#'
#' @examples data(ScotchWhiskey)
#' lapply(ScotchWhiskey,ncol)
#' ScotchWhiskey$nbChar
#' ScotchWhiskey$listW  ## attr(ScotchWhiskey$listW,"class")
#' names(ScotchWhiskey)
#' names(ScotchWhiskey$dist)
#' ##
#' plotWhiskey <- function(main) {
#'     plot(x=ScotchWhiskey$geo@coords[,1L]/1000,
#'          xlab="Eastings (km)",
#'          y=ScotchWhiskey$geo@coords[,2L]/1000,
#'          ylab="Northings (km)",
#'          main=main,
#'          type="n",asp=1)
#'     apply(
#'         ScotchWhiskey$neighbor@data,1L,
#'         function(X,coords) {
#'             segments(
#'                 coords[X[1L],1L]/1000,
#'                 coords[X[1L],2L]/1000,
#'                 coords[X[2L],1L]/1000,
#'                 coords[X[2L],2L]/1000
#'             )
#'         },
#'         coords=ScotchWhiskey$geo@coords
#'     )
#'     invisible(NULL)
#' }
#' ##
#' plotWhiskey("Scotch whiskey: peat nose")
#' cols <- c("blue","orange")
#' points(ScotchWhiskey$geo@coords/1000,pch=21L,
#'        bg=cols[ScotchWhiskey$nose[,"peat"]+1L])
#' legend(x=50,y=1000,legend=c("Has a peat nose","Has no peat nose"),
#'        pch=21L,pt.bg=rev(cols))
#' ##
#' plotWhiskey("Scotch whiskey: soft body")
#' cols <- c("red","green")
#' points(ScotchWhiskey$geo@coords/1000,pch=21L,
#'        bg=cols[ScotchWhiskey$body[,"soft"]+1L])
#' legend(x=50,y=1000,legend=c("Has a soft body","Has no soft body"),
#'        pch=21L,pt.bg=rev(cols))
#' ##
#' plotWhiskey("Scotch whiskey: spicy palate")
#' cols <- c("red","green")
#' points(ScotchWhiskey$geo@coords/1000,pch=21L,
#'        bg=cols[ScotchWhiskey$palate[,"spice"]+1L])
#' legend(x=50,y=1000,legend=c("Has a spicy palate","Has no spicy palate"),
#'        pch=21L,pt.bg=rev(cols))
#' ##
#' plotWhiskey("Scotch whiskey: sweet finish")
#' cols <- c("red","green")
#' points(ScotchWhiskey$geo@coords/1000,pch=21L,
#'        bg=cols[ScotchWhiskey$finish[,"sweet"]+1L])
#' legend(x=50,y=1000,legend=c("Has a sweet finish","Has no sweet finish"),
#'        pch=21L,pt.bg=rev(cols))
#' ##
#' ### To visualize (part of) the distance matrices:
#' as.matrix(ScotchWhiskey$dist$nose)[1:5,1:5]
#' as.matrix(ScotchWhiskey$dist$body)[1:5,1:5]
#' as.matrix(ScotchWhiskey$dist$palate)[1:5,1:5]
#' as.matrix(ScotchWhiskey$dist$finish)[1:5,1:5]
#' ##
#' ### The data tables:
#' ScotchWhiskey$colour
#' head(ScotchWhiskey$nose)
#' head(ScotchWhiskey$body)
#' head(ScotchWhiskey$palate)
#' head(ScotchWhiskey$finish)
NULL
