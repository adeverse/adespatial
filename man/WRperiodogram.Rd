\name{WRperiodogram}
\alias{WRperiodogram}
\alias{plot.WRperio}
\title{ Whittaker-Robinson periodogram }

\description{
  Whittaker-Robinson periodogram for univariate series of quantitative data.
}

\usage{
WRperiodogram(x, T1 = 2, T2, nperm = 499, nopermute, mult = c("sidak",
"bonferroni"), print.time = FALSE)
\method{plot}{WRperio}(x, prog = 1, alpha = 0.05, line.col = "red", ...)
}

\arguments{
  \item{x}{ A vector of quantitative values, with class \code{numeric}, for function \code{WRperiodogram}, or an output object of WRperiodogram for function \code{plot}. }
  \item{T1}{ First period included in the calculation (default: \code{T1=2}). }
  \item{T2}{ Last period included in the calculation (default: \code{T2=n/2}). }
  \item{nperm}{ Number of permutations for the tests of significance.  }
  \item{nopermute}{ List of item numbers that should not be permuted; see Details (default: no items should be excluded from the permutations). }
  \item{mult}{ Correction method for multiple testing. Choices are "bonferroni" and "sidak" (default: \code{mult="bonferroni"}). }
  \item{print.time}{ Print the computation time. Useful when planning the analysis of a long data series with a high number of permutations. Default: \code{print.time=FALSE}. }
  \item{prog}{ \code{prog=1} (default): use the original p-values in the plot.
	\code{prog=2}: use the p-values corrected for multiple testing.
	\code{prog=3}: progressive correction of the multiple tests. }
  \item{alpha}{ Significance level for the plot; p-values smaller than or equal to alpha are represented by black symbols. Default: \code{alpha=0.05}. }
  \item{line.col}{ Colour of the lines between symbols in the graph (default: \code{line.col="red"}). }
  \item{\dots}{ Other graphical arguments passed to this function. }
} 

\details{ The Whittaker-Robinson periodogram (Whittaker and Robinson, 1924) identifies periodic components in a vector of quantitative data. The data series must contain equally-spaced observations (i.e. constant lag) along a transect in space or through time. The vector may contain missing observations, represented by NA, in reasonable amount, e.g. up to a few percent of the total number of observations. The periodogram statistic used in this function is the standard deviation of the means of the columns of the Buys-Ballot table (Enright, 1965). The method is also described in Legendre & Legendre (2012, Section 12.4.1). Missing values (NA) are handled by skipping the NA values when computing the column means of the Buys-Ballot table.

The data must be stationary before computation of the periodogram. Stationarity is violated when there is a trend in the data or when they were obtained under contrasting environmental or experimental conditions. Users should at least test for the presence of a significant linear trend in the data (using linear regression); if a significant trend is identified, it can be removed by computing regression residuals.

The \code{nopermute} option allows users to include a list of items numbers that should not be permuted, whether the observations are NA or zero values. This option should not be used in routine work. It is intended for special situations where observations could not be made at some points along the space or time series because that was impossible. For example, in a spatial data series along a river, if points fall on emerging rocks or on islands, no observation of phytoplankton could have been made at those points. For the permutation test, values at these positions (NA or 0) should not be permuted with values at points where observations were possible. 

The graph produced by the \code{plot} function shows the periodogram statistics and their significance following a permutation test, with periods in the abscissa. The p-values may be corrected for multiple testing using either the Bonferroni or the Sidak correction, which can be applied to all values in the correlogram uniformly, or following a progressive correction.

A progressive correction means that for the first periodogram statistic, the p-value is tested against the alpha significance level without any correction; for the second statistic, the p-value is corrected for 2 simultaneous tests; and so forth until the k-th statistic, where the p-value is corrected for k simultaneous tests. This approach solves the problem of "where to stop interpreting a periodogram"; one goes on as long as significant values emerge, considering the fact that the tests become progressively more conservative.

In the Whittaker-Robinson periodogram, harmonics of a basic period are often found to be also significant.

The permutation tests, which can take a bit of time in very large jobs, can be interrupted by issuing an \code{Escape} command. One can also click the \code{STOP} button at the top of the R console.
}

\value{ 
The function produces an object of class \code{WRperio} containing a table with the following columns: 

\item{Period}{ period number; }
\item{WR.stat}{	periodogram statistic; } 
\item{p-value}{	p-value after permutation test; }
\item{p.corrected}{	p-value corrected for multiple tests, using the Bonferroni or Sidak method; }
\item{p.corr.prog}{ p-value after progressive correction. }

When the p-values cannot be computed because of a very high proportion of missing values in the data, values of 99 are posted in the last three columns of the output table.
}

\author{ Pierre Legendre and Guillaume Guenard, Universite de Montreal}

\references{

  Enright, J. T. 1965. The search for rhythmicity in biological time-series. 
  Journal of Theoretical Biology 8: 426-468.

  Legendre, P. and L. Legendre. 2012. Numerical ecology, 3rd English
  edition. Elsevier Science BV, Amsterdam.

  Sarrazin, J., D. Cuvelier, L. Peton, P. Legendre and P. M. Sarradin. 2014. High-resolution dynamics of a deep-sea hydrothermal mussel assemblage monitored by the EMSO-Açores MoMAR observatory. Deep-Sea Research I 90: 62-75. (Recent application in oceanography)

  Whittaker, E. T. and G. Robinson. 1924. The calculus of observations – 
  A treatise on numerical mathematics. Blackie & Son, London.
}

\examples{ 
###
### 1. Numerical example of Subsection 12.4.1 of Legendre and Legendre (2012)

test.vec <- c(2,2,4,7,10,5,2,5,8,4,1,2,5,9,6,3)

# Periodogram with permutation tests of significance
res <- WRperiodogram(test.vec)
plot(res)

### 2. Simulated data

periodic.component <- function(x,T,c) cos((2*pi/T)*(x+c))

n <- 500   # corresponding to 125 days, 4 observations per day
# Generate a lunar cycle, 29.5 days (T=118)
moon <- periodic.component(1:n, 118, 59)
# Generate a circadian cycle (T=4)
daily <- periodic.component(1:n, 4, 0)
# Generate a tidal cycle (T=2)
tide <- periodic.component(1:n, 2, 0)

# Periodogram of the lunar component only 
res.moon <- WRperiodogram(moon, nperm=0)
res.moon <- WRperiodogram(moon, T2=130, nperm=99)
par(mfrow=c(1,2))
plot(moon)
plot(res.moon, prog=1)

# Add the three components, plus a random normal error term
var <- 5*moon + daily + tide + rnorm(n, 0, 0.5)
# Draw a graph of a portion of the data series
par(mfrow=c(1,2))
plot(var[1:150], pch=".", cex=1)
lines(var[1:150])

# Periodogram of 'var'
res.var <- WRperiodogram(var, T2=130, nperm=99)
plot(res.var, prog=1, line.col="blue")
# Find position of the maximum value of this periodogram
which(res.var[,2] == max(res.var[,2]))

# Replace 10% of the 500 data by NA
select <- sort(sample(1:500)[1:50])
var.na <- var
var.na[select] <- NA
res.var.na <- WRperiodogram(var.na, T2=130, nperm=99)
plot(res.var.na, prog=1)

### 3. Data used in the examples of the documentation file of function afc() of {stats}
# Data file "ldeaths"; time series, 6 years x 12 months of deaths in UK hospitals
ld.res.perio <- WRperiodogram(ldeaths, nperm=499)
par(mfrow=c(1,2))
plot(ld.res.perio, prog=1) # Graph with no correction for multiple testing
plot(ld.res.perio, prog=3) # Graph with progressive correction
acf(ldeaths)   # acf() results, for comparison
}

\keyword{ multivariate }
