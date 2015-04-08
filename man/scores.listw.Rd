\name{scores.listw}
\alias{scores.listw}
\alias{mem}
\alias{orthobasis.listw}
\title{Function to compute Moran's eigenvectors of a listw object}
\description{
  This functions compute eigenvector's of a doubly centered spatial
  weighting matrix. Corresponding eigenvalues are linearly related to
  Moran's index of spatial autocorrelation.
}
\usage{
scores.listw(listw, wt = rep(1, length(listw$neighbours)), MEM.autocor =
c("non-null", "all","positive", "negative"))
orthobasis.listw(listw, wt = rep(1, length(listw$neighbours)), MEM.autocor =
c("non-null", "all","positive", "negative"))
mem(listw, wt = rep(1, length(listw$neighbours)), MEM.autocor =
c("non-null", "all","positive", "negative"))
}

\arguments{
  \item{listw}{ An object of the class \code{listw} created by functions
    of the \code{spdep} package
  }
  \item{wt}{A vector of weights. It is used to orthogonalize the
  eigenvectors. It could be useful if MEM are used in weighted
  regression or canonical correspondence analysis
}
\item{MEM.autocor}{ A string indicating if all MEMs must be returned or
  only those corresponding to non-null, positive or negative
  autocorrelation
}
}
\details{
  Testing the nullity of eigenvalues is based on E(i)/E(1) where E(i) is
  i-th eigenvalue and E(1) is the maximum absolute value of eigenvalues
  
}
\value{
  An object of class \code{orthobasis}
}
\references{
  Dray, S., Legendre, P., and Peres-Neto, P. R. (2006). Spatial modeling: a
  comprehensive framework for principal coordinate analysis of neighbor
  matrices (PCNM).\emph{Ecological Modelling} \bold{196}, 483--493.

  
  Griffith D. A. (1996) Spatial autocorrelation and eigenfunctions of
  the geographic weights matrix accompanying geo-referenced
  data. \emph{Canadian Geographer} \bold{40}, 351--367.
}
\author{Stéphane Dray}



\seealso{ \code{\link[spdep]{nb2listw}}}
\examples{
if(require("ade4", quietly = TRUE) & require("spdep", quietly = TRUE)){
data(oribatid)
nbtri <- tri2nb(as.matrix(oribatid$xy))
sc.tri <- scores.listw(nb2listw(nbtri, style = "B"))
summary(sc.tri)
}
if(require("adegraphics", quietly = TRUE)){
s.value(oribatid$xy,sc.tri[,1:9])
}
}

\keyword{ spatial }
