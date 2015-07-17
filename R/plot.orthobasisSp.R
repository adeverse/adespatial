#' Function to display Moran's Eigenvector Maps (MEM) and other spatial 
#' orthogonal bases
#' 
#' This function allow to plot or map orthogonal bases
#' @author Stéphane Dray \email{stephane.dray@@univ-lyon1.fr}
#' @param x an object of class \code{orthobasisSp}
#' @param SpORcoords either a \code{Spatial*} object or a \code{matrix} with 
#'   geographic coordinates
#' @param pos an integer indicating the position of the environment where the 
#'   data are stored, relative to the environment where the function is called. 
#'   Useful only if ‘storeData’ is ‘FALSE’
#' @param storeData a logical indicating if the data should be stored in the 
#'   returned object. If ‘FALSE’, only the names of the data arguments are 
#'   stored
#' @param plot a logical indicating if the graphics is displayed
#' @param match.ID a logical indicating if names of geographic entities match 
#'   rownames of the \code{orthobasisSp} object
#' @param \dots additional graphical parameters (see ‘adegpar’ and 
#'   ‘trellis.par.get’)
#'   
#' @importFrom adegraphics s.Spatial
#'   
#' @examples
#' if(require("ade4", quietly = TRUE) & require("spdep", quietly = TRUE)){
#' data(mafragh)
#' me <- mem(nb2listw(mafragh$nb))
#' 
#' if(require("adegraphics", quietly = TRUE)){
#' plot(me[,1:6], mafragh$xy)
#' plot(me[,1:6], mafragh$Spatial) 
#' }
#' }
#'         
#' @export
plot.orthobasisSp <- function(x, SpORcoords, pos = -1, storeData = TRUE, plot = TRUE, match.ID = FALSE, ...){
    if(missing(SpORcoords)){
        NextMethod()
    } else {
        if(length(grep("Spatial", class(SpORcoords))) > 0){
            SpObjectCall <- substitute(SpORcoords)
        } else {
            SpObjectCall <- call("SpatialPoints", substitute(as.matrix(SpORcoords)))
            SpORcoords <- do.call("SpatialPoints", list(substitute(as.matrix(SpORcoords))))           
        }
        x <- substitute(data.frame(x))
        if(inherits(SpORcoords, what = "SpatialPoints"))
            SpObjectCall <- call("SpatialPointsDataFrame", coords = SpObjectCall, data = x, match.ID = match.ID)
        if(inherits(SpORcoords, what = "SpatialPolygons"))
            SpObjectCall <- call("SpatialPolygonsDataFrame", Sr = SpObjectCall, data = x, match.ID = match.ID)
        if(inherits(SpORcoords, what = "SpatialGrid"))
            SpObjectCall <- call("SpatialGridDataFrame", grid = SpObjectCall, data = x)

        object <- do.call("s.Spatial", list (spObj = SpObjectCall, plot = FALSE, storeData = TRUE, pos = pos - 2, ...))
        if(plot)
            print(object)
        invisible(object)
    }        
}
