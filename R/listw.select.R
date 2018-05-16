#'Function to optimize the selection of a spatial weighting matrix and select
#'the best subset of eigenvectors
#'
#'\code{listw.select} computes MEM variables (i.e., eigenvectors of a doubly
#'centered spatial weighting matrix) for various definitions of spatial
#'weighting matrices (SWM) and optimizes the selection of the SWM and of a
#'subset of MEM variables. The optimization is done by maximizing the adjusted
#'R-squared (R2) or by minimizing the residual spatial autocorrelation. The
#'function controls the type I error rate by accounting for the number of tests
#'performed. This function combine calls to the functions \code{scores.listw} and
#'\code{mem.select}. The list of candidate SWMs
#'can easily be generated using \code{\link{listw.candidates}}. 
#'
#'
#'@details While the selection of the SWM is the most critical step of the
#'  spatial eigenvector-based methods (Dray et al. 2006), Bauman et al. (2018a)
#'  showed that optimizing the choice of the SWM led to inflated type I error
#'  rates if an explicit control of the number of SWMs tested was not applied.
#'  The function \code{listw.select} therefore applies a Sidak correction (Sidak
#'  1967) for multiple tests to the p-value of the global test of each SWM
#'  (i.e., the model integrating the whole set of spatial predictors). The Sidak
#'  correction is computed as: \eqn{P_corrected = 1 - (1 - P)^n}, where \eqn{n}
#'  is the number of tests performed, \eqn{P} is the observed p-value, and
#'  \eqn{P_corrected} is the new p-value after the correction. The p-value is
#'  first computed using \code{nperm} permutations and then corrected according to the total number
#'  of SWMs tested (if \code{p.adjust = TRUE}). Although the function can be
#'  run without this correction, using the default
#'  value is strongly recommended to avoid inflated
#'  type I error rates (Bauman et al. 2018a).
#'
#'  As a consequence of the p-value correction, the
#'  significance threshold decreases as the number of SWMs increases, hence
#'  leading to a trade-off between the gain of accuracy and the power loss. 
#'
#'  The optimization criterion of the SWM performed by \code{listw.select} is
#'  either based on the maximization of the adjusted R2 of all the generated
#'  spatial eigenvectors (also referred to as spatial predictors or MEM
#'  variables) (\code{method = "global"}), or is based on an optimized subset of
#'  eigenvectors (\code{method= "FWD"} or \code{"MIR"}). If the objective is
#'  only to optimize the selection of the SWM, without the intervention of the
#'  selection of a subset of predictors within each SWM (\code{method=
#'  "global"}), then the best SWM will be the one maximizing the adjusted global
#'  R2, that is, the R2 of the model of \code{x} against the whole set of
#'  generated MEM variables. If \code{MEM.autocor = "all"}, then the adjusted R2 is
#'  computed separately on the MEM associated to positive and negative
#'  eigenvalues (hereafter positive and negative MEM variables, respectively),
#'  and the SWM yielding the highest sum of the the two R2 values is selected.
#'  Using \code{method = "global"} may be interesting when the complete set of MEM
#'  variables will be used, like in Moran spectral randomizations (Wagner and
#'  Dray 2015) or when using smoothed MEM (Munoz 2009).
#'
#'  If the MEM variables will be further used in a model including actual
#'  predictors (e.g. environmental), then a subset of spatial eigenvectors will
#'  need to be selected before proceeding to further analyses to avoid model
#'  overfitting and a loss of statistical power to detect the contribution of
#'  the environment to the variability of the response data (Griffith 2003, Dray
#'  et al. 2006, Blanchet et al. 2008, Peres-Neto and Legendre 2010, Diniz-Filho
#'  et al. 2012). In this case, the optimization of the SWM should be based on
#'  the best subset of spatial predictors within each candidate SWM. Although
#'  several eigenvector selection approaches have been proposed to select this
#'  best subset of eigenvectors, Bauman et al. (2018b) showed that two main
#'  procedures should be preferred, depending on the underlying objective.
#'
#'  The most powerful and accurate selection method, in terms of R2 estimation,
#'  is the forward selection with double stopping criterion of Blanchet et al.
#'  (2008). This method (\code{method = "FWD"}, default option) should be
#'  preferred when the objective is to capture as accurately as possible the
#'  spatial patterns of \code{x}. If the objective is to optimize the detection
#'  of the spatial patterns in the residuals of a model of the response
#'  variable(s) against a set of environmental predictors, for instance, then
#'  \code{x} should be the model residuals, and \code{FWD = "forward"}. This
#'  can be interesting if the objective is to optimize the detection of residual
#'  spatial patterns once the effect of the environmental predictors has been
#'  removed.
#'
#'  If however there is no interest in exploring accurate spatial structures in
#'  the residuals of the model of \code{x} against the environmental predictors,
#'  and the objective is only to remove the spatial autocorrelation from the
#'  model residuals with a small number of spatial predictors, then accuracy is
#'  not as important and one should focus mainly on the number of spatial
#'  predictors (Bauman et al. 2018b). In this case, using \code{method = "MIR"} is
#'  the most adapted choice. This optimization criterion selects the smallest
#'  subset of spatial predictors necessary to capture all the spatial
#'  autocorrelation in \code{x} (see function \code{\link{mem.select}}. It has
#'  the advantage to maintain the standard errors of the actual predictor
#'  coefficients as low as possible.
#'
#'  If \code{method = "FWD"}, \code{listw.select} performs the forward
#'  selection on the significant SWMs and selects among these the SWM for which
#'  the forward-selected subset of spatial eigenvectors yields the highest
#'  adjusted R2. If \code{method = "MIR"}, \code{listw.select} performs the MIR
#'  selection on all the significant candidate SWMs, and selects the best SWM as
#'  the one with the smallest number of MIR-selected spatial eigenvectors. If
#'  two or more SWMs present the same smallest number of predictors, then the
#'  selection is made among them on the basis of the residual Moran's I. If \code{MEM.autocor = "all"}, the
#'  optimization criteria described above are applied on the sum of the adjusted
#'  R2 or number of selected spatial eigenvectors, for \code{method = "FWD"}
#'  and \code{"MIR"}, respectively.
#'
#'  Note that the MIR criterion of optimization can only be used for a
#'  univariate \code{x}, as the Moran's I is a univariate index. If \code{x} is
#'  multivariate, then the best criterion is the forward selection (see Bauman
#'  et al. 2018b).
#'
#'  \code{mem.select} works exactly like \code{listw.select} but for a single
#'  SWM. Its purpose is to test the SWM and, if significant, to optimize the
#'  selection of a subset of spatial predictors within it.
#'
#'@param x Vector, matrix, or dataframe of the response variable(s)
#'@param candidates A list of SWMs of the class \code{listw};
#'  \code{candidates} can be created by \code{listw.candidates}
#' @param MEM.autocor Sign of the spatial eigenvectors to generate; "positive",
#'   "negative", or "all", for positively, negatively autocorrelated
#'   eigenvectors, or both, respectively; default is "positive"
#' @param method Criterion to select the best subset of MEM variables. Either
#'   \code{forward} (default option), \code{"MIR"} (for univariate \code{x}
#'   only), or \code{"global"} (see \code{Details})
#' @param MEM.all A logical indicating if the complete set of MEM variables for the best model
#'   should be returned
#' @param nperm Number of permutations to perform the tests; Default is 999
#' @param alpha Significance threshold value for the tests; Default is 0.05
#' @param p.adjust A logical indicating wheter the p-value of the global test performed on each SWM
#'  should be corrected for multiple tests (TRUE) or not (FALSE); default is
#'  \code{TRUE}
#'
#'@return \code{listw.select} returns a list that contains:
#'  \describe{ \item{candidates}{A data.frame that summarizes the results on all SWMs}
#'  \item{best.id}{The index and name of the best SWM} \item{best}{The results 
#'  for the best SWM as returned by \code{mem.select}} }
#'
#'
#'@author Bauman David (\email{dbauman@@ulb.ac.be} or
#'  \email{davbauman@@gmail.com}) and Stéphane Dray
#'
#'@seealso \code{\link{listw.candidates}}, \code{\link{mem.select}},
#'  \code{\link{scores.listw}}
#'
#'@references Bauman D., Fortin M-J, Drouet T. and Dray S. (2018a) Optimizing
#'  the choice of a spatial weighting matrix in eigenvector-based methods.
#'  Ecology
#'
#'  Bauman D., Drouet T., Dray S. and Vleminckx J. (2018b) Disentangling good
#'  from bad practices in the selection of spatial or phylogenetic eigenvectors.
#'  Ecography, 41, 1--12
#'
#'  Blanchet G., Legendre P. and Borcard D. (2008) Forward selection of
#'  explanatory variables. Ecology, 89(9), 2623--2632
#'
#'  Diniz-Filho J.A.F., Bini L.M., Rangel T.F., Morales-Castilla I. et al.
#'  (2012) On the selection of phylogenetic eigenvectors for ecological
#'  analyses. Ecography, 35, 239--249
#'
#'  Dray S., Legendre P. and Peres-Neto P. R. (2006) Spatial modeling: a
#'  comprehensive framework for principal coordinate analysis of neighbor
#'  matrices (PCNM). Ecological Modelling, 196, 483--493
#'
#'  Griffith D. (2003) Spatial autocorrelation and spatial filtering: gaining
#'  understanding through theory and scientific visualization. Springer, Berlin
#'
#'  Munoz, F. 2009. Distance-based eigenvector maps (DBEM) to analyse
#'  metapopulation structure with irregular sampling. Ecological Modelling, 220,
#'  2683–-2689
#'
#'  Peres-Neto P. and Legendre P. (2010) Estimating and controlling for spatial
#'  structure in the study of ecological communities. Global Ecology and
#'  Biogeography, 19, 174--184
#'
#'  Sidak Z. (1967) Rectangular confidence regions for the means of multivariate
#'  normal distributions. Journal of the American Statistical Association,
#'  62(318), 626--633
#'
#'  Wagner H., Dray S. (2015). Generating spatially constrained null models for
#'  irregularly spaced data using Moran spectral randomization methods. Methods
#'  in Ecology and Evolution, 6, 1169--1178
#'
#'@keywords spatial
#'
#' @examples
#' if(require(spdep)) {
#' ### Create a grid of 15 x 15:
#' grid <- expand.grid(x = seq(1, 15, 1), y = seq(1, 15, 1))
#' ### Generate a response variable Y structured at broad scale by linear combination of
#' ### the first three MEM variables to which a normal noise is added:
#' nb <- cell2nb(nrow = 15, ncol = 15, "queen")
#' nb2 <- nb2listw(nb, style = "B")
#' MEM <- scores.listw(nb2, MEM.autocor = "positive")
#' # Degree of spatial autocorrelation:
#' intensity <- 0.8
#' Y_space <- scale(MEM[, 1] + MEM[, 2] + MEM[, 3]) * intensity
#' Y_noise <- scale(rnorm(n = nrow(MEM), mean = 0, sd = 1)) * (1 - intensity)
#' Y <- Y_space + Y_noise
#' ### Y is sampled in 100 randomly-chosen sites of the grid:
#' sample <- sample(c(1:nrow(grid)), 100, replace = FALSE)
#' xy <- grid[sample, ]
#' Y_sampled <- Y[sample]
#' ### The function listw.candidates is used to build the spatial weighting matrices that
#' ### we want to test and compare (with the listw.select function). We test a Gabriel's graph,
#' ### a minimum spanning tree, and a distance-based connectivity defined by a threshold
#' ### distance corresponding to the smallest distance keeping all sites connected (i.e.,
#' ### the defaut value of d2; see help of function listw.candidates).
#' ### These connectivity matrices are then either not weighted (binary weighting), or
#' ### weighted by the linearly decreasing function (see help of the function listw.candidates):
#' candidates <- listw.candidates(coord = xy, nb = c("gab", "mst"), weights = c("binary", "flin"))
#' ### Number of candidate W matrices generated:
#' nbw <- length(candidates)
#' ### Significance threshold value after p-value correction (Sidak correction):
#' 1 - (1 - 0.05)^(1/nbw)
#' ### Optimization of the selection of the SWM among the candidates generated above,
#' ### using the corrected significance threshold calculated above for the global tests:
#' W_sel <- listw.select(Y_sampled, candidates, MEM.autocor = "positive", method = "FWD",
#'                     p.adjust = TRUE)
#' ### Some characteristics of the best spatial model:
#' # Best SWM:
#' W_sel$best.id
#' # Selected subset of spatial predictor within the best SWM:
#' W_sel$best$MEM.select
#' nrow(W_sel$best$summary)
#' # Corrected p-value of the global test of the best SWM:
#' W_sel$best$global.test$Pvalue
#' # Adjusted R2 of the subset of spatial predictors selected within the chosen SWM:
#' max(W_sel$best$summary$R2Adj)
#' # p-values of all the tested W matrices:
#' W_sel$candidates$Pvalue
#' # Adjusted R2 of the subset of spatial predictors selected for all the significant
#' # W matrices:
#' W_sel$candidates$R2Adj.select
#'
#' # See Appendix S3 of Bauman et al. 2018a for more extensive examples and illustrations.
#' }
#'
#'@export listw.select
#'@importFrom stats na.omit

"listw.select" <- function(x, 
    candidates, 
    MEM.autocor = c("positive", "negative", "all"), 
    method = c("FWD", "MIR", "global"),
    MEM.all = FALSE, nperm = 999, alpha = 0.05, 
    p.adjust = TRUE) {
    
    
    method <- match.arg(method)
    MEM.autocor <- match.arg(MEM.autocor)
    
    ntest <- ifelse(p.adjust, length(candidates), 1)
    res.tmp <- lapply(candidates, mem.select, x = x, MEM.autocor = MEM.autocor, 
        method = method, MEM.all = MEM.all, nperm = nperm, alpha = alpha, ntest = ntest)
    
    res <- data.frame(row.names = names(candidates))
    
    .getstat <- function(list.res, MEM.autocor, method){
        ## get the statistic of the global test
        df <- data.frame(row.names = names(list.res))
        prefix <- ifelse(method == "MIR", "I", "R2Adj")
        if(MEM.autocor == "all")
            suffix <- c(".pos", ".neg") 
        else 
            suffix <- NULL
        if(MEM.autocor == "all"){
            df[,1] <- sapply(list.res, function(x) x$global.test$positive$obs)
            df[,2] <- sapply(list.res, function(x) x$global.test$negative$obs)
        } else {
            df[,1] <- sapply(list.res, function(x) x$global.test$obs)
        }
        names(df) <- paste(prefix, suffix, sep = "")
        return(df) 
    }
    
    
    
    .getpvalue <- function(list.res, MEM.autocor){
        ## get the pvalue of the global test
        df <- data.frame(row.names = names(list.res))
        if(MEM.autocor == "all")
            suffix <- c(".pos", ".neg") 
        else 
            suffix <- NULL
        if(MEM.autocor == "all"){
            df[,1] <- sapply(list.res, function(x) x$global.test$positive$pvalue)
            df[,2] <- sapply(list.res, function(x) x$global.test$negative$pvalue)
        } else {
            df[,1] <- sapply(list.res, function(x) x$global.test$pvalue)
        }
        names(df) <- paste("Pvalue", suffix, sep = "")
        return(df) 
    }
    
    .getinfoselect <- function(list.res, method){
        ## get the information of the best subset of MEM
        df <- data.frame(row.names = names(list.res))
        df[,1] <- sapply(list.res, function(x) ifelse(is.null(x$summary), NA, nrow(x$summary)))
        if(method == "MIR"){
            df[,2] <- sapply(list.res, function(x) ifelse(is.null(x$summary), NA, x$summary[nrow(x$summary), "Iresid"]))
            names(df)<- c("N.var", "I.select")
        } else {
            df[,2] <- sapply(list.res, function(x) ifelse(is.null(x$summary), NA, x$summary[nrow(x$summary), "AdjR2Cum"]))
            names(df)<- c("N.var", "R2Adj.select")
        }
        return(df) 
    }
    
    
    res <- cbind(res, .getstat(res.tmp, MEM.autocor = MEM.autocor, method = method))
    res <- cbind(res, .getpvalue(res.tmp,  MEM.autocor = MEM.autocor))
    res <- cbind(res, .getinfoselect(res.tmp, method = method))
    
    ## Select the best model among candidates
    if (length(res$N.var) > sum(is.na(res$N.var))) {
        if(method == "FWD") {
            best <- which.max(res$R2Adj.select) 
        } else if(method == "global") {
            best <- which.max(ifelse(MEM.autocor == "all", res$R2Adj.pos + res$R2Adj.neg, res$R2Adj))
        } else if(method == "MIR") {
            best <- which(res$N.var == min(na.omit(res$N.var)))
            if(length(best) > 1)
                ## if different models have the same nummber of variables
                ## minimizes  residual autocorrelation
                best <- best[which.min(abs(res$I.select[best]))] 
        }
        
        best.id <- best
        names(best.id) <- names(res.tmp)[best]
        return(list(candidates = res, best.id = best.id, best = res.tmp[[best]]))   
    } else {
        return(list(candidates = res))
        
    }
}
