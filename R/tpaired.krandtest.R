#' Paired t-tests of differences between T1 and T2 for each species
#' 
#' This function computes paired t-tests for each species, for abundances observed 
#' at time 1 (T1) and time 2 (T2). The test is one-tailed in the direction of the sign 
#' (+ or -) of the t statistic. 
#' 
#' @param mat1 site-by-species data at time T1 (data.frame or matrix).
#' @param mat2 site-by-species data at time T2 (data.frame or matrix).
#' @param nperm  Number of permutations. Use 999, 9999, or more, 
#'    to allow for correction of p-values for multiple tests.
#' @param list.all If \code{FALSE}, the output matrix \code{$t.tests} only lists t.test  
#'   results for species with \code{t.stat} not 0;
#'   If \code{TRUE}, the output matrix \code{$t.tests} lists t.test results for all 
#'   species; when \code{t.stat} is 0, the p-values in the output table (\code{p.param} 
#'   and \code{p.perm}) receive codes -999; \code{Sign(T1-T2)} receives the value 0.
#' 
#' @details The species that do not vary in either data set are discarded before 
#' calculation of the paired t-tests begins.
#' 
#' p-values should be corrected for multiple testing. Use function \code{p.adjust} of 
#' \code{stats}: p.adjust(res$t.test$p.param)  or  p.adjust(res$t.test$p.perm)
#' Correction methods "holm" (default) and "hochberg" are fine for this type of analysis.
#' 
#' @return \itemize{
#' \item A table with species in rows and 6 columns: 
#' "mean(T1-T2)","t.stat","p.param","p.perm","p<=0.05","Sign(T1-T2)"
#' The parametric and permutational p-values are not corrected for multiple tests.
#' A star is shown in column "p<=0.05" if the parametric p-value is <= 0.05.
#' \item A list of names of the species tested; their t statistics were not 0.
#' \item A list of names of the species not tested because their t-statistics 
#' were 0. }
#' 
#' @author Pierre Legendre \email{pierre.legendre@@umontreal.ca}
#' 
#' @references Legendre, P. 2019. A temporal beta-diversity index to identify sites that 
#' have changed in exceptional ways in space-time surveys. \emph{Ecology and Evolution} 
#' (in press).
#' 
#' van den Brink, P. J. & C. J. F. ter Braak. 1999. Principal response curves: analysis of 
#' time-dependent multivariate responses of biological community to stress. 
#' \emph{Environmental Toxicology and Chemistry} 18: 138-148.
#' 
#' @seealso \code{tpaired.randtest}
#' 
#' @examples
#' 
#' if(require("vegan", quietly = TRUE)) {
#' 
#' ## Invertebrate communities subjected to insecticide treatment.
#' 
#' ## As an example in their paper on Principal Response Curves (PRC), van den Brink & ter 
#' ## Braak (1999) used observations on the abundances of 178 invertebrate species 
#' ## (macroinvertebrates and zooplankton) subjected to treatments in 12 mesocosms by the 
#' ## insecticide chlorpyrifos. The mesocosms were sampled at 11 occasions. The data, 
#' ## available in the {vegan} package, are log-transformed species abundances, 
#' ## y.tranformed = loge(10*y+1).
#' 
#' ## The data of survey #4 will be compared to those of survey #11 in this example.  
#' ## Survey #4 was carried out one week after the insecticide treatment, whereas the 
#' ## fauna of the mesocosms was considered to have fully recovered from the treatment 
#' ## at the time of survey #11.
#' 
#' data(pyrifos)
#' 
#' ## The mesocosms had originally been attributed at random to the treatments. However,  
#' ## to facilitate presentation of the results, they will be listed here in order of 
#' ## increased insecticide doses: {0, 0, 0, 0, 0.1, 0.1, 0.9, 0.9, 6, 6, 44, 44} 
#' ## micro g/L.
#' 
#' survey4.order = c(38,39,41,47,37,44,40,46,43,48,42,45)
#' 
#' survey11.order = c(122,123,125,131,121,128,124,130,127,132,126,129)
#' 
#' ## Paired t-tests of differences between survey.4 and survey.11 for the p species
#' 
#' res <- tpaired.krandtest(pyrifos[survey4.order,],pyrifos[survey11.order,])
#' 
#' }
#' 
#' @export 

tpaired.krandtest <- function(mat1, mat2, nperm=99, list.all=FALSE) {
    epsilon <- .Machine$double.eps
    n <- nrow(mat1)
    p <- ncol(mat1)
    if (nrow(mat2) != n) stop("Unequal number of rows")
    if (ncol(mat2) != p) stop("Unequal number of species")
    
    # Check species names in files mat1 and mat2
    sp.diff <- which(colnames(mat1) != colnames(mat2))
    if (length(sp.diff) > 0) 
        cat("Warning: the following species names differ between mat1 & mat2:\n",sp.diff,"\n\n")
    
    # Select the species that have variances > epsilon. Discard the other species.
    comm.gr <- rbind(mat1, mat2)
    tmp <- apply(comm.gr, 2, var)
    sel.sp <- which(tmp > epsilon)
    pp <- length(sel.sp)
    cat((p - pp),"species were eliminated because they did not vary in the combined data set\n")
    cat(pp,"species retained:\n")
    mat1.sel <- mat1[,sel.sp]
    mat2.sel <- mat2[,sel.sp]
    
    res <- matrix(0,pp,4)
    sp.names <- rownames(res) <- colnames(mat1.sel)
    
    Star <- NA
    Keep <- NA
    Reject <- NA
    for (j in 1:pp) {
        vec1 <- mat1.sel[,j]
        vec2 <- mat2.sel[,j]
        # if(any(abs(vec1-vec2) > 0)) Keep = c(Keep,j) else Reject = c(Reject,j)
        cat(j," ")      # Print ID of species being processed
        tmp <- t.test(vec1, vec2, paired = TRUE)
        if (tmp$estimate == 0) {
            Reject <- c(Reject,j) 
            res[j, 3:4] <- c(-999,-999)
            Star <- c(Star, " ")
        } else {
            Keep = c(Keep,j) 
            if (sign(tmp$estimate) < 0)
                tail <- "less"
            else if (sign(tmp$estimate) >= 0) 
                tail <- "greater"
            t.res <- tpaired.randtest(vec1, vec2, nperm = nperm, alternative = tail, silent = TRUE)
            res[j,1:4] <- c(t.res$estim, t.res$t.ref, t.res$p.param, t.res$p.perm)
            Star <- c(Star, ifelse(t.res$p.perm > 0.05, " ","*"))
        }
    }
    cat("\n\n")
    Keep <- Keep[-1]
    Reject <- Reject[-1]
    tmp <- length(Reject)
    if (tmp > 0) 
        cat(tmp,"species not tested because t.stat = 0. See 'No_test' output list\n\n")
    res <- data.frame(res,Star[-1],sign(res[,1]))
    colnames(res) <- c("mean(T1-T2)","t.stat","p.param","p.perm","p<=0.05","Sign(T1-T2)")
    if (list.all) 
        out <- list(t.tests = res, Tested = sp.names[Keep], No_test = sp.names[Reject])
    else 
        out <- list(t.tests = res[Keep,], Tested = sp.names[Keep], No_test = sp.names[Reject])
    out
}
