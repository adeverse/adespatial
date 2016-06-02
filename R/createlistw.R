#' Interactive tool to generate R code that creates a spatial weighting matrix
#' 
#' @export createlistw
#' @importFrom shiny runApp
#' @author Stéphane Dray \email{stephane.dray@@univ-lyon1.fr}
#' @seealso \code{\link{chooseCN}}
#' @examples
#' \dontrun{
#' ## a matrix or an object of class 'Spatial*' should be in the global environment
#' xy <- matrix(rnorm(50), 25)
#' createlistw()
#' }
createlistw <- function (){
    runApp(system.file('createlistw', package = 'adespatial'))
}