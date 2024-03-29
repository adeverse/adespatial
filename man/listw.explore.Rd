% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/listw.explore.R
\name{listw.explore}
\alias{listw.explore}
\title{Interactive tool to generate R code that creates a spatial weighting matrix}
\usage{
listw.explore()
}
\value{
No return value
}
\description{
Interactive tool to generate R code that creates a spatial weighting matrix
}
\examples{
if(interactive()){
## a matrix or an object of class 'Spatial*' should be in the global environment
xy <- matrix(rnorm(50), 25)
listw.explore()
}
}
\seealso{
\code{\link{chooseCN}}
}
\author{
Stéphane Dray \email{stephane.dray@univ-lyon1.fr}
}
