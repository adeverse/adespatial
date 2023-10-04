#' Interactive tool to generate R code that creates a spatial weighting matrix
#' 
#' @export listw.explore
#' @importFrom shiny runApp
#' @author Stéphane Dray \email{stephane.dray@@univ-lyon1.fr}
#' @seealso \code{\link{chooseCN}}
#' @returns No return value
#' @examples
#' if(interactive()){
#' ## a matrix or an object of class 'Spatial*' should be in the global environment
#' xy <- matrix(rnorm(50), 25)
#' listw.explore()
#' }
listw.explore <- function (){
    runApp(system.file('listw.explore', package = 'adespatial'))
}