#' Permutational paired t-test
#'
#' This function computes a permutation test of comparison of the means of two
#' paired vectors (related samples). For each object, repeated permutations are
#' restricted to the two related observations.
#'
#' @param vec1,vec2 The two data vectors to be compared.
#' @param nperm Number of permutations. Use 999, 9999, or more, in real-case applications.
#' @param alternative c("two.sided", "less", "greater"). Default value:
#'   "two.sided".
#' @param silent If \code{FALSE}, calculation results are printed to the R
#'   console. If \code{TRUE} calculation results are not printed to R console
#'   (e.g. for simulations, or for use inside another function).
#'
#' @return A list containing the following results: \itemize{ \item
#'   \code{estim}: mean of the differences \item \code{t.ref}: reference value
#'   of the t-statistic \item \code{p.param}: parametric p-value \item
#'   \code{p.perm}: permutational p-value \item \code{nperm}: number of
#'   permutations }
#'
#' @author Pierre Legendre \email{pierre.legendre@@umontreal.ca}. Permutation
#'   code improved by Guillaume Blanchet \email{guillaume.blanchet@@usherbrooke.ca}.
#'
#' @references Zar, J. H. 1999. \emph{Biostatistical analysis. 4th edition.}
#'   Prentice Hall, New Jersey.
#'
#' @seealso \code{\link[stats]{t.test}}
#'
#' @examples
#'
#' ## Deer leg length, data from Zar (1999, p. 162).
#'
#' deer <- matrix(c(142,140,144,144,142,146,149,150,142,148,138,136,147,139,143,141,143,
#' 145,136,146),10,2)
#'
#' rownames(deer) <- paste("Deer",1:10,sep=".")
#'
#' colnames(deer) <- c('Hind.leg', 'Fore.leg')
#'
#' res <- tpaired.randtest(deer[,1], deer[,2])   # Two-tailed test by default
#'
#' ## Compare the results to:  res2 = t.test(deer[,1], deer[,2], paired=TRUE)
#'
#' @importFrom stats t.test rbinom
#' @export tpaired.randtest
#'   

tpaired.randtest <- function(vec1, vec2, nperm = 99, alternative = "two.sided", silent = FALSE) {
    n1 <- length(vec1)
    n2 <- length(vec2)
    if (n1 != n2) stop("The two vectors have different lengths. They cannot be paired.")
    
    alt <- match.arg(alternative, c("two.sided", "less", "greater"))
    
    res <- t.test(vec1, vec2, paired = TRUE, alternative = alt)
    t.ref <-  res$statistic
    
    # Print these first results
    if (!silent) { 
        cat('\nt-test comparing the means of two related samples','\n','\n')
        cat('Number of objects:',n1,'\n')
        cat('Mean of the differences:',res$estimate,'\n')
        cat('t statistic (paired observations):',t.ref,'\n')
        cat('95 percent confidence interval of t:',res$conf.int,'\n')
        cat('Degrees of freedom:',res$parameter,'\n')
        cat('Alternative hypothesis:',alt,'\n')
        cat('Prob (parametric):',res$p.value,'\n')
    }
    
    # Perform the permutation test
    # Permutations are restricted to the two related observations for each object.
    nPGE <- 1
    
    for (i in 1:nperm) {
        mat <- cbind(vec1,vec2)
        topermute <- rbinom(n1,1,0.5)
        mat[topermute == 1,] <- mat[topermute == 1, 2:1]
        
        res.perm <- t.test(mat[,1], mat[,2], paired = TRUE, alternative = alt)
        t.perm <- res.perm$statistic
        
        if (alt == "two.sided") if (abs(t.perm) >= abs(t.ref) ) nPGE <- nPGE + 1
        if (alt == "less")      if (t.perm <= t.ref) nPGE <- nPGE + 1
        if (alt == "greater")   if (t.perm >= t.ref) nPGE <- nPGE + 1
    }
    
    # Compute and print the permutational p-value
    P <- nPGE / (nperm + 1)
    if (!silent) cat('Prob (',nperm,'permutations):', formatC(P, digits = 5, width = 7,format = "f"),'\n')
    #
    return(list(estim = res$est[[1]], t.ref = t.ref, p.param = res$p.value, p.perm = P, nperm = nperm))
}
