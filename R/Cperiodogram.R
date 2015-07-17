#'Contingency periodogram
#'
#'Function to compute a contingency periodogram for a univariate series of 
#'qualitative data
#'
#'The contingency periodogram of Legendre et al. (1981) identifies periodic 
#'components in qualitative data vectors. The vector may contain classes of a 
#'qualitative variable or the classes obtained by hierarchical clustering or 
#'partitioning of a multivariate data table. The method is also described in 
#'Legendre & Legendre (2012). The optional graph produced by the function shows 
#'the following information:
#'
#'\itemize{ \item In red: the B statistics (information in common). \item In 
#'blue: Confidence limits for B without correction. \item In green: 
#'Bonferroni-corrected limits of the confidence intervals. \item In black: 
#'Confidence limits with progressive Bonferroni correction. }
#'
#'@param x a qualitative variable (\code{factor})
#'@param T1 first period included in the calculations (default: T1 = 2)
#'@param T2 last period included in the calculations (default: T2 = n/2)
#'@param nperm Number of permutations for the chi-square test. For chi-square 
#'  tests using the chi- square distribution, use the default \code{nperm=NULL}
#'@param alpha significance level for computation of the confidence limits
#'@param graph a logical indicating if a graph is requested, by default 
#'  \code{TRUE}.
#'  
#'@author Pierre Legendre \email{pierre.legendre@@umontreal.ca}
#'  
#'@return A table with the statistics for the selected periods: \itemize{ \item
#'  Wilks’ chi-square statistic (Wilks.chisq) \item information in common (B), 
#'  \item degrees of freedom (df), \item p-value (prob) } Confidence interval
#'  limits: \itemize{ \item critical value of B without correction (B.crit), 
#'  \item critical value of B with Bonferroni correction based on the number of
#'  periods studied in the periodogram (B.crit.Bonf), \item critical value of B
#'  with progressive Bonferroni correction (B.prog.Bonf). }
#'  
#'  @references
#'  Legendre, L., M. Fréchette & P. Legendre. 1981. The contingency periodogram: a method of
#' identifying rhythms in series on nonmetric ecological data. Journal of Ecology 69: 965-979.
#'
#' Legendre, P. and Legendre, L. 2012. Numerical Ecology. 3rd English ed. Elsevier, Amsterdam
#' 
#' @examples
#' # Data from the numerical example of Subsection 12.4.2 of Legendre and Legendre (2012).
#' test.vec <- c(1,1,2,3,3,2,1,2,3,2,1,1,2,3,3,1)
#' # Periodogram with tests using the chi-square distribution
#' res <- Cperiodogram(test.vec)
#' # Periodogram with permutation tests
#' res <- Cperiodogram(test.vec, nperm=2000, graph=FALSE)
#' 
#' @export
Cperiodogram <- function(x, T1 = 2, T2 = NULL, nperm = NULL, alpha = 0.05, graph = TRUE)
# Contingency periodogram (Legendre et al. 1981).
#
# License: GPL-2 
# Author:: Pierre Legendre, November 2010
{
### Internal functions
###
Wilks.chisq <- function(Obs, Exp)
	{
	rem <- which(Obs==0)
	if(length(rem)==0)
		{
		res <- 2*sum(Obs * log(Obs/Exp))
		} else {
		res <- 2*sum(Obs[-rem] * log(Obs[-rem]/Exp[-rem]))
		}
	}

plo <- function(out)
	{
	plot(out[,1], out[,3], col="red", xlab="Period", 
		ylab="Information in common (B)", ylim=c(0, max(out[,3], out[,7])))
	lines(out[,1], out[,3], col="red")     # B statistics
	# points(out[,1], out[,6], col="blue")
	lines(out[,1], out[,6], col="blue")    # Ordinary confidence limits
	# points(out[,1], out[,7], col="green")
	lines(out[,1], out[,7], col="green")   # Bonferroni-corrected conf. limits
	# points(out[,1], out[,8], col="black")
	lines(out[,1], out[,8], col="black")   # Progressive Bonferroni conf. limits
	}
###
### End internal functions

if(!is.factor(x)) if(length(which(x <= 0)) > 0) 
	stop("x contains 0 or negative levels")
x <- as.factor(x)
nlev.x <- nlevels(x)
n <- length(x)
T.sup <- n %/% 2
if(is.null(T2)) T2 <- T.sup
if(T2 > T.sup) T2 <- T.sup
nk <- T2-T1+1
if(is.null(nperm)) {
	simulate <- FALSE
	} else {
	simulate <- TRUE
	}
out <- matrix(NA,nk,8)
colnames(out) <-c("Period","Wilks.chisq","B","df","prob","B.crit","B.crit.Bonf","B.prog.Bonf")
for(k in T1:T2) {
	# cat("Period =",k,"\n")
	vec <- as.factor(rep(1:k, ceiling(n/k)))
	tab <- table(x,vec[1:n])
	chisq.out <- chisq.test(tab, simulate.p.value=simulate, B=nperm)
	chisq <- Wilks.chisq(chisq.out$observed, chisq.out$expected)
	B <- chisq/(2*n)
	if(simulate) {df <- (nrow(tab)-1)*(ncol(tab)-1)} 
		else {df <- chisq.out$parameter}
	out[(k-T1+1),1] <- k
	out[(k-T1+1),2] <- chisq
	out[(k-T1+1),3] <- B
	out[(k-T1+1),4] <- df
	out[(k-T1+1),5] <- chisq.out$p.value
	out[(k-T1+1),6] <- qchisq(alpha, df, lower.tail=FALSE)/(2*n)
	out[(k-T1+1),7] <- qchisq(alpha/nk, df, lower.tail=FALSE)/(2*n)
	out[(k-T1+1),8] <- qchisq(alpha/(k-T1+1), df, lower.tail=FALSE)/(2*n)
	}
if(graph) plo(out)
out
}