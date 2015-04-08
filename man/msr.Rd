\name{msr}
\alias{msr}
\title{
Moran spectral randomization
}
\description{
This function allows to generate spatially-constrained random variables
preserving the global autocorrelation (Moran's I) and the spatial
structures at multiple scales. Multiscale property is defined by the
power spectrum (i.e. decomposition of the variance of the original
variables) on a basis of orthonormal eigenvectors (Moran's Eigenvector
Maps, MEM). The function provides methods for univariate randomization,
joint randomization of a group of variables while keeping within-group
correlations fixed and univariate randomization with a fixed correlation
between original data and randomized replicates.
}
\usage{
msr(x, listwORorthobasis, nrepet = 99, method = c("pair", "triplet",
"singleton"), cor.fixed, nmax = 100, simplify = TRUE)
}

\arguments{
  \item{x}{a \code{vector}, a \code{matrix} or a \code{data.frame} with
  the original variables. If \code{NCOL(x) > 1}, then the joint
  randomization procedure that preserves the correlations among
  variables is used.
}
  \item{listwORorthobasis}{an object of the class \code{listw} (spatial
  weights) created by the functions of the \pkg{spdep} package or an
  object of class \code{orthobasis}
}
  \item{nrepet}{an \code{integer} indicating the number of replicates}
  \item{method}{an character specifying which algorithm should be used
  to produce spatial replicates (see Details).
}
  \item{cor.fixed}{if not missing, the level of correlation between the
  original variable and its randomized replicates
}
  \item{nmax}{the number of trials used in the "triplet" procedure.
}
  \item{simplify}{A logical value. If \code{TRUE}, the outputs for
  univariate procedures are returned in a matrix where each column
  corresponds to a replicate. If \code{FALSE}n a \code{list} is returned.
}
}
\details{
Three procedures are implemented in the function. The "pair" procedure
  is the more general as it can be applied in the three cases
  (univariate, univariate with fixed correlation and multivariate). This
  procedure preserves the power spectrum by pair of MEMs but not
  strictly the global autocorrelation level (Moran's I). The "singleton"
  procedure can be used for univariate and multivariate cases. It
  preserves strictly the global level of autocorrelation and the power
  spectrum. The "triplet" procedure can only be applied in the
  univariate case. It preserves the power  spectrum by triplet of MEMs
  and strictly the global autocorrelation level.
  
}
\value{
Either a matrix (if \code{simplify} is
\code{TRUE}) or a list with randomized replicates. 
}
\references{
Wagner, H.H. and Dray S. (2015) Generating spatially-constrained null models for irregularly spaced data using Moran spectral randomization methods
}
\author{
Stephane Dray \email{stephane.dray@univ-lyon1.fr} and Helene H Wagner \email{helene.wagner@utoronto.ca}
}




\seealso{
\code{\link{scores.listw}}, \code{\link[spdep]{nb2listw}}
}
\examples{
library(spdep)
x1 <- matrix(rnorm(81*5), nrow = 81)
lw1 <- nb2listw(cell2nb(9, 9))

moran.mc(x1[,1], lw1, 2)$statistic

## singleton
x1.1 <- msr(x1[,1], lw1, nrepet = 9, method = "singleton")
apply(x1.1, 2, function(x) moran.mc(x, listw = lw1, nsim = 2)$statistic)

## triplet
x1.2 <- msr(x1[,1], lw1, nrepet = 9, method = "triplet")
apply(x1.2, 2, function(x) moran.mc(x, listw = lw1, nsim = 2)$statistic)

## pair
x1.3 <- msr(x1[,1], lw1, nrepet = 9, method = "pair")
apply(x1.3, 2, function(x) moran.mc(x, listw = lw1, nsim = 2)$statistic)

## pair with cor.fixed
x1.4 <- msr(x1[,1], lw1, nrepet = 9, cor.fixed = 0.5)
apply(x1.4, 2, function(x) moran.mc(x, listw = lw1, nsim = 2)$statistic)
cor(x1[,1], x1.4)

## pair preserving correlations for multivariate data
x1.5 <- msr(x1, lw1, nrepet = 9, cor.fixed = 0.5)
cor(x1)
lapply(x1.5, cor)

apply(x1, 2, function(x) moran.mc(x, listw = lw1, nsim = 2)$statistic)
apply(x1.5[[1]], 2, function(x) moran.mc(x, listw = lw1, nsim = 2)$statistic)

## singleton preserving correlations for multivariate data
x1.6 <- msr(x1, lw1, nrepet = 9, method = "singleton")
cor(x1)
lapply(x1.6, cor)

apply(x1, 2, function(x) moran.mc(x, listw = lw1, nsim = 2)$statistic)
apply(x1.6[[1]], 2, function(x) moran.mc(x, listw = lw1, nsim = 2)$statistic)
}

\keyword{ spatial }

