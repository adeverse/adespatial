% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tpaired.krandtest.R
\name{tpaired.krandtest}
\alias{tpaired.krandtest}
\title{Paired t-tests of differences between T1 and T2 for each species}
\usage{
tpaired.krandtest(mat1, mat2, nperm = 99, list.all = FALSE, opt = 2)
}
\arguments{
\item{mat1}{site-by-species data at time T1 (data.frame or matrix).}

\item{mat2}{site-by-species data at time T2 (data.frame or matrix).}

\item{nperm}{Number of permutations. Use 999, 9999, or more, 
to allow for correction of p-values for multiple tests. In this documentation 
file, \code{nperm} was set to 99 only for demonstration of the function.}

\item{list.all}{If \code{FALSE}, the output matrix \code{$t.tests} only lists t.test
results for species with \code{t.stat} not 0;
If \code{TRUE}, the output matrix \code{$t.tests} lists t.test results for all 
species; when \code{t.stat} is 0, the p-values in the output table (\code{p.param} 
and \code{p.perm}) receive codes -999; \code{Sign(T1-T2)} receives the value 0.}

\item{opt}{Compute the difference (T1-T2) if \code{opt=1}, or the difference (T2-T1) 
if \code{opt=2} (default). \code{opt=2} is the default because this produces the 
easiest-to-read output for the comparison of successive surveys in time series analysis 
of community composition data. \code{opt=1} is also available to users with different 
needs. The only differences between the two sets of results are the signs of the 
'\code{mean(T1-T2)} and \code{t.stat} statistics. The original version of this function 
only offered opt=1. opt=2 was added in April 2025 in answer to a user's request.}
}
\value{
\itemize{
\item A table with species in rows and 6 columns: 
"mean(T1-T2)" or "mean(T2-T1)" ,"t.stat","p.param","p.perm","p<=0.05","Sign(T1-T2)" 
or "Sign(T2-T1)".
The parametric and permutational p-values are not corrected for multiple tests.
This correction can be applied to a column of p-values in the function's output file 
(see previous paragraph).
A star is shown in column \code{p<=0.05} if the parametric p-value is <= 0.05. 
\item A list of names of the species tested; their t statistics were not 0.
\item A list of names of the species not tested because their t-statistics were 0. }
}
\description{
This function computes paired t-tests for each species, for abundances observed 
at time 1 (T1) and time 2 (T2). The test is one-tailed in the direction of the sign 
(+ or -) of the t statistic.
}
\details{
The species that do not vary in either data set are discarded before 
calculation of the paired t-tests begins.

p-values should be corrected for multiple testing. Use function \code{p.adjust} of 
\code{stats}: \code{p.adjust(res$t.test$p.param)} or \code{p.adjust(res$t.test$p.perm)}.
Correction methods \code{"holm"} (default) and \code{"hochberg"} are both fine for this 
type of analysis.
}
\examples{

if(require("vegan", quietly = TRUE)) {

## Invertebrate communities subjected to insecticide treatment.

## As an example in their paper on Principal Response Curves (PRC), van den Brink & ter 
## Braak (1999) used observations on the abundances of 178 invertebrate species 
## (macroinvertebrates and zooplankton) subjected to treatments in 12 mesocosms by the 
## insecticide chlorpyrifos. The mesocosms were sampled at 11 occasions. The data, 
## available in the {vegan} package, are log-transformed species abundances, 
## y.tranformed = loge(10*y+1).

## The data of survey #4 will be compared to those of survey #11 in this example.  
## Survey #4 was carried out one week after the insecticide treatment, whereas the 
## fauna of the mesocosms was considered to have fully recovered from the treatment 
## at the time of survey #11.

data(pyrifos)

## The mesocosms had originally been attributed at random to the treatments. However,  
## to facilitate presentation of the results, they will be listed here in order of 
## increased insecticide doses: {0, 0, 0, 0, 0.1, 0.1, 0.9, 0.9, 6, 6, 44, 44} 
## micro g/L.

survey4.order = c(38,39,41,47,37,44,40,46,43,48,42,45)

survey11.order = c(122,123,125,131,121,128,124,130,127,132,126,129)

## Paired t-tests of differences between survey.4 and survey.11 for the p species

res <- tpaired.krandtest(pyrifos[survey4.order,],pyrifos[survey11.order,], opt=2)

## Smaller problem where only species 15 to 19 are analysed

mat1 = pyrifos[survey4.order,]
mat2 = pyrifos[survey11.order,]
res2 <- tpaired.krandtest(mat1[,15:19], mat2[,15:19], opt=2)

\dontrun{  # More permutations to obtain a more precise estimate of the p-value

res2 <- tpaired.krandtest(mat1[,15:19], mat2[,15:19], nperm=9999, opt=2)
}

}

}
\references{
Legendre, P. 2019. A temporal beta-diversity index to identify sites that 
have changed in exceptional ways in space-time surveys. \emph{Ecology and Evolution} 
(in press).

van den Brink, P. J. & C. J. F. ter Braak. 1999. Principal response curves: analysis of 
time-dependent multivariate responses of biological community to stress. 
\emph{Environmental Toxicology and Chemistry} 18: 138-148.
}
\seealso{
\code{tpaired.randtest}
}
\author{
Pierre Legendre \email{pierre.legendre@umontreal.ca}
}
