#' Function to optimize the selection of a spatial weighting matrix and select the best 
#' subset of eigenvectors
#' 
#' \code{listw.select} computes MEM variables (i.e., eigenvectors of a doubly centered spatial 
#' weighting matrix) for various definitions of spatial weighting matrices (SWM) and optimizes 
#' the selection of the SWM and of a subset of MEM variables. 
#' The optimization is done by maximizing the adjusted R-squared (R2) or by minimizing the 
#' residual spatial autocorrelation. The function controls the type I error rate by accounting 
#' for the number of tests performed. 
#' Function \code{mem.select} computes MEM variables for a single definition of SWM, test its 
#' significances, and optimizes the selection of MEM variables within the SWM using the same
#' criteria than \code{listw.select}.
#' These functions combine calls to the functions \code{scores.listw} and \code{forward.sel} 
#' or \code{\link{MEM.moransel}}.
#' The list of candidate SWMs can easily be generated using \code{\link{listw.candidates}}.
#' The significance of each SWM is tested by 9999 permutations with the function 
#' \code{anova.cca} (package \code{vegan}).
#' 
#' @details While the selection of the SWM is the most critical step of the spatial
#' eigenvector-based methods (Dray et al. 2006), Bauman et al. (2018a) showed that 
#' optimizing the choice of the SWM led to inflated type I error rates if an
#' explicit control of the number of SWMs tested was not applied. The function
#' \code{listw.select} therefore applies a Sidak correction (Sidak 1967) for multiple tests to
#' the p-value of the global test of each SWM (i.e., the model integrating the 
#' whole set of spatial predictors). The Sidak correction is computed as: 
#'\eqn{P_corrected = 1 - (1 - P)^n}, where \eqn{n} is the number of tests performed, \eqn{P}
#' is the observed p-value, and \eqn{P_corrected} is the new p-value after the correction.
#' The p-value is first computed by 9999 permutations with the function \code{anova.cca} 
#' (package \code{vegan}), and is then corrected according to the total number of SWMs 
#' tested (if \code{correction = TRUE}). Although the function can be run without this
#' correction (\code{correction = FALSE}), using the default value (\code{correction = TRUE})
#' is strongly recommended to avoid inflated type I error rates (Bauman et al. 2018a).
#' 
#' As a consequence of the necessity of the p-value correction, the significance threshold
#' decreases as the number of SWMs increases, hence leading to a trade-off between the gain 
#' of accuracy and the power loss. See \code{Details} of the function 
#' \code{\link{listw.candidates}} for recommendations about the selection of a 
#' suitable set of candidate SWMs. 
#' 
#' The optimization criterion of the SWM performed by \code{listw.select} is either based 
#' on the maximization of the adjusted R2 of all the generated spatial eigenvectors (also 
#' referred to as spatial predictors or MEM variables) (\code{crit = "global"}), or is based 
#' on an optimized subset of eigenvectors (\code{crit = "forward"} or \code{"MIR"}).
#' If the objective is only to optimize the selection of the SWM, without the intervention 
#' of the selection of a subset of predictors within each SWM (\code{crit = "global"}), then 
#' the best SWM will be the one maximizing the adjusted global R2, that is, the R2 of the model
#' of \code{x} against the whole set of generated MEM variables. If \code{autocor = "all"}, 
#' then the adjusted R2 is computed separately on the MEM associated to positive and negative 
#' eigenvalues (hereafter positive and negative MEM variables, respectively), and the SWM 
#' yielding the highest sum of the the two R2 values is selected.
#' Using \code{crit = "global"} may be interesting when the complete set of MEM variables will
#' be used, like in Moran spectral randomizations (Wagner and Dray 2015) or when using 
#' smoothed MEM (Munoz 2009). 
#' 
#' If the MEM variables will be further used in a model including actual predictors 
#' (e.g. environmental), then a subset of spatial eigenvectors will need to
#' be selected before proceeding to further analyses to avoid model overfitting and a loss of 
#' statistical power to detect the contribution of the environment to the 
#' variability of the response data (Griffith 2003, Dray et al. 2006, Blanchet et al. 2008,
#' Peres-Neto and Legendre 2010, Diniz-Filho et al. 2012). In this case, the optimization
#' of the SWM should be based on the best subset of spatial predictors within each
#' candidate SWM. Although several eigenvector selection approaches have been proposed to 
#' select this best subset of eigenvectors, Bauman et al. (2018b) showed that two main 
#' procedures should be preferred, depending on the underlying objective.
#' 
#' The most powerful and accurate selection method, in terms of R2 estimation, is the 
#' forward selection with double stopping criterion of Blanchet et al. (2008). This method
#' (\code{crit = "forward"}, default option) should be preferred when the objective is to
#' capture as accurately as possible the spatial patterns of \code{x}.
#' If the objective is to optimize the detection of the spatial patterns in the residuals of
#' a model of the response variable(s) against a set of environmental predictors, for instance,
#' then \code{x} should be the model residuals, and \code{crit = "forward"}. This can be 
#' interesting if the objective is to optimize the detection of residual spatial patterns
#' once the effect of the environmental predictors has been removed.
#' 
#' If however there is no interest in exploring accurate spatial structures in the residuals
#' of the model of \code{x} against the environmental predictors, and the objective is only
#' to remove the spatial autocorrelation from the model residuals with a small number of
#' spatial predictors, then accuracy is not as important and one should focus mainly on the
#' number of spatial predictors (Bauman et al. 2018b).
#' In this case, using \code{crit = "MIR"} is the most adapted choice. This optimization 
#' criterion selects the smallest subset of spatial predictors necessary to capture all the 
#' spatial autocorrelation in \code{x} (see function \code{\link{MEM.moransel}}. It has the 
#' advantage to maintain the standard errors of the actual predictor coefficients as low as 
#' possible. 
#' 
#' If \code{crit = "forward"}, \code{listw.select} performs the forward selection on the 
#' significant SWMs and selects among these the SWM for which the forward-selected
#' subset of spatial eigenvectors yields the highest adjusted R2. 
#' If \code{crit = "MIR"}, \code{listw.select} performs the MIR selection on all the 
#' significant candidate SWMs, and selects the best SWM as the one with the smallest number 
#' of MIR-selected spatial eigenvectors. If two or more SWMs present the same smallest number 
#' of predictors, then the selection is made among them on the basis of the adjusted R2 of the 
#' selected eigenvectors (as for a same number of spatial eigenvectors, a higher R2 indicates 
#' a more accurate description of the spatial structure in the residuals (Bauman et al. 2018a)).
#' If \code{autocor = "all"}, the optimization criteria described above are applied on the
#' sum of the adjusted R2 or number of selected spatial eigenvectors, for 
#' \code{crit = "forward"} and \code{"MIR"}, respectively.
#' 
#' Note that the MIR criterion of optimization can only be used for a univariate \code{x}, 
#' as the Moran's I is a univariate index. If \code{x} is multivariate, then the best 
#' criterion is the forward selection (see Bauman et al. 2018b).
#' 
#' \code{mem.select} works exactly like \code{listw.select} but for a single SWM. Its purpose
#' is to test the SWM and, if significant, to optimize the selection of a subset of spatial
#' predictors within it.
#'
#' @aliases listw.select mem.select  
#' @param x Vector, matrix, or dataframe of the response variable(s)
#' @param candidates A list of one or more SWMs of the class \code{listw}; \code{candidates} 
#' can be a list of several SWM (when using \code{listw.select}), but it can also be a list 
#' of one SWM or a single \code{listw} object (when using \code{mem.select})
#' @param autocor Sign of the spatial eigenvectors to consider; "positive", "negative", 
#' or "all", for positively, negatively autocorrelated eigenvectors, or both, respectively; 
#' default is "positive"
#' @param crit Criterion to select the best SWM among the set of significant candidate 
#' matrices. Either \code{forward} (default option), \code{"MIR"} (for univariate \code{x} 
#' only), or \code{"global"} (see \code{Details})
#' @param alpha Uncorrected predefined significance value; default is 0.05
#' @param correction wheter the p-value of the global test performed on each SWM should
#' be corrected for multiple tests (TRUE) or not (FALSE); default is \code{TRUE}
#' @param summary Logical; Whether a message summarizing the results should be returned
#' or not; The R2 value returned in the message depends on the argument \code{crit} (R2 of 
#' the subset of selected predictors or R2 of the global model); Default is \code{TRUE}
#' 
#' @return \code{listw.select} returns two lists containing several features of the SWM selected
#' by the optimization procedure, and general features of all the SWMs generated and tested: 
#' 
#' The first list, \code{$all}, contains: 
#' 
#' \describe{
#' \item{MEM.all}{A list of all the complete set of MEM variables (as generated by 
#' \code{scores.listw}) for the candidate SWMs provided to the function, in the same order as 
#' the SWMs in \code{candidates}.}
#' \item{pval}{The p-values (corrected by the Sidak correction if \code{correction = TRUE})
#' of the global models for each SWM. If \code{autocor = "all"}, two values per SWM (positive 
#' and negative MEM models).}
#' \item{R2.global}{The adjusted R2 of the complete set of MEM variables of each SWM.
#' R2.global is provided only for the significant SWMs. If \code{autocor = "all"}, two values 
#' per SWM (positive If \code{autocor = "all"}, two values per SWM (positive 
#' and negative MEM models).and negative MEM models).}
#' }
#' 
#' The second list, \code{$best}, contains the features of the best model: 
#' 
#' \describe{
#' \item{MEM.all}{Dataframe of the complete set of eigenvectors of the best SWM.}
#' \item{MEM.select}{Dataframe of the subset of eigenvectors selected by the criterion
#' \code{crit}. This set of spatial predictors is the one to be used in further analyses 
#' (e.g. variation partitioning, ordinary least squares, GLM) (see \code{Details} section). 
#' If \code{crit = "global"}, then \code{MEM.select} contains the entire set of MEM variables 
#' generated for the corresponding SWM. In this case, \code{MEM.select} and \code{MEM.all}
#' are the same.} 
#' \item{listw}{Element of \code{candidates} corresponding to the best SWM. It is an 
#' object of class \code{listw}.}
#' \item{MEM.AdjR2Cum}{Only for \code{crit = "forward"}; Vector of the cumulative adjusted 
#' R2 of the eigenvectors selected by the forward selection (\code{MEM.select}).}
#' \item{name}{Name of the best SWM ("Connectivity matrix_Weighting matrix") if the 
#' list of candidate SWMs was generated by the function \code{\link{listw.candidates}}. 
#' Otherwise, corresponds to the index of the best SWM in the list \code{candidates} 
#' provided to \code{listw.select}.} 
#' \item{pval}{Global p-value of the best SWM (after multiple-test correction, if 
#' the \code{correction} argument = TRUE). If \code{autocor = "all"}, two values (positive and 
#' negative MEM models).}
#' \item{R2.global}{The adjusted R2 of the complete set of MEM variables of the best W 
#' matrix. If \code{autocor = "all"}, two values (positive and negative MEM models).}
#' \item{R2.select}{Adjusted R2 of \code{MEM.select}. If \code{crit = "global"}, then 
#' \code{R2.select} is equal to \code{R2.global}. If \code{autocor = "all"}, two values 
#' (positive and negative MEM models).}
#' \item{NbVar}{Number of spatial predictors selected by the forward selection. If 
#' \code{autocor = "all"}, two values (positive and negative MEM models).}
#' \item{bestw_index}{Index of the best SWM in the list \code{candidates} provided 
#' to \code{listw.select}.}
#' }
#' 
#' The function \code{mem.select} only returns one list, corresponding to the \code{$best} list
#' of \code{listw.selec}, but without the objects \code{name} and \code{bestW_index}.
#' 
#' @author Bauman David \email{dbauman@@ulb.ac.be} or \email{davbauman@@gmail.com}
#' 
#' @seealso \code{\link{listw.candidates}}, \code{\link{MEM.moransel}}, 
#' \code{\link{scores.listw}}, \code{\link[vegan]{varpart}}
#' 
#' @references Bauman D., Fortin M-J, Drouet T. and Dray S. (2018a) Optimizing the choice of 
#' a spatial weighting matrix in eigenvector-based methods. Ecology
#' 
#' Bauman D., Drouet T., Dray S. and Vleminckx J. (2018b) Disentangling good from bad 
#' practices in the selection of spatial or phylogenetic eigenvectors. Ecography, 41, 1--12
#' 
#' Blanchet G., Legendre P. and Borcard D. (2008) Forward selection of explanatory variables.
#' Ecology, 89(9), 2623--2632
#' 
#' Diniz-Filho J.A.F., Bini L.M., Rangel T.F., Morales-Castilla I. et al. (2012) On the 
#' selection of phylogenetic eigenvectors for ecological analyses. Ecography, 35, 239--249
#' 
#' Dray S., Legendre P. and Peres-Neto P. R. (2006) Spatial modeling: a comprehensive 
#' framework for principal coordinate analysis of neighbor matrices (PCNM). Ecological 
#' Modelling, 196, 483--493
#' 
#' Griffith D. (2003) Spatial autocorrelation and spatial filtering: gaining understanding 
#' through theory and scientific visualization. Springer, Berlin
#' 
#' Munoz, F. 2009. Distance-based eigenvector maps (DBEM) to analyse metapopulation structure 
#' with irregular sampling. Ecological Modelling, 220, 2683–-2689
#' 
#' Peres-Neto P. and Legendre P. (2010) Estimating and controlling for spatial structure 
#' in the study of ecological communities. Global Ecology and Biogeography, 19, 174--184
#' 
#' Sidak Z. (1967) Rectangular confidence regions for the means of 
#' multivariate normal distributions. Journal of the American Statistical Association, 62(318),
#' 626--633
#' 
#' Wagner H., Dray S. (2015). Generating spatially constrained null models for irregularly 
#' spaced data using Moran spectral randomization methods. Methods in Ecology and Evolution, 
#' 6, 1169--1178
#' 
#' @keywords spatial
#' 
#' @examples
#' if(require(spdep)) {
#' ### Create a grid of 15 x 15:
#' grid <- expand.grid(x = seq(1, 15, 1), y = seq(1, 15, 1))
#' ### Generation of a response variable Y structured at broad scale by linear combination of
#' ### the first three MEM variables to which a normal noise is added:
#' nb <- cell2nb(nrow = 15, ncol = 15, "queen")
#' nb2 <- nb2listw(nb, style = "B")
#' MEM <- scores.listw(nb2, MEM.autocor = "positive")
#' # Degree of spatial autocorrelation:
#' intensity <- 0.4
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
#' candidates <- listw.candidates(coord = xy, gab = TRUE, mst = TRUE, binary = TRUE, 
#'                                flin = TRUE)
#' ### Number of candidate W matrices generated:
#' nbw <- length(candidates)
#' ### Significance threshold value after p-value correction (Sidak correction):
#' 1 - (1 - 0.05)^(1/nbw)
#' ### Optimization of the selection of the SWM among the candidates generated above, 
#' ### using the corrected significance threshold calculated above for the global tests:
#' W_sel <- listw.select(Y_sampled, candidates, autocor = "positive", crit = "forward", 
#'                     correction = TRUE)
#' ### Some characteristics of the best spatial model:
#' # Best SWM:
#' W_sel$best$name
#' # Selected subset of spatial predictor within the best SWM:
#' W_sel$best$MEM.select
#' W_sel$best$NbVar
#' # Corrected p-value of the global test of the best SWM: 
#' W_sel$best$pval
#' # Adjusted R2 of the subset of spatial predictors selected within the chosen SWM:
#' W_sel$best$R2.select
#' # p-values of all the tested W matrices:
#' W_sel$all$pval
#' # Adjusted R2 of the subset of spatial predictors selected for all the significant
#' # W matrices:
#' W_sel$all$R2.select
#' 
#' # See Appendix S3 of Bauman et al. 2018a for more extensive examples and illustrations.
#' }
#' 
#' @export listw.select mem.select 
#' @importFrom vegan rda anova.cca RsquareAdj
#' @importFrom stats na.omit

"listw.select" <- function(x, 
                           candidates, 
                           autocor = c("positive", "negative", "all"), 
                           crit = c("forward", "MIR", "global"),
                           alpha = 0.05, 
                           correction = TRUE,
                           summary = TRUE) {
  
  x <- as.data.frame(x)
  if (any(is.na(x))) stop ("NA entries in x")
  
  crit <- match.arg(crit)
  if (crit == "MIR") 
    if (ncol(x) != 1) stop ("MIR criterion can only be used for a univariate x")
  
  if (correction == TRUE) sidak <- "corrected" else sidak <- "uncorrected"
  autocor <- match.arg(autocor)

  # **********************************************************************************  
  MEM.test <- function(response = x, 
                       mem, 
                       selection = crit,
                       c = autocor, 
                       d = nbtest, 
                       thresh = alpha,
                       corr = correction) {
    # MEM.test is an internal function that tests the significance of a SWM while taking
    # into consideration the total number of W matrices tested in listw.select (if corr = TRUE).
    # This total number of tests is used to apply a correction to the p-value in order not to
    # inflate the type I error rate (if corr = TRUE). If the tested SWM is significant 
    # at the corrected threshold value of significance, then either a model selection is 
    # performed using Blanchet et al.'s forward selection (2008) (crit = "forward") or the
    # MIR selection (crit = "MIR"), or only the adjusted R2 of the global model is computed
    # (crit = "global").
    pval <- anova.cca(rda(response, mem), permutations = 9999)$Pr[1]
    if (correction == TRUE) pval <- 1-(1-pval)^d    # Sidak correction 
    if (c == "all") pval <- 1-(1-pval)^2            # Sidak correction (autocor= "all") 
    if (pval <= thresh) {  
      R2 <- RsquareAdj(rda(response, mem))$adj.r.squared
      if (selection == "forward") {
        class <- class(try(fsel <- forward.sel(response, mem, adjR2thresh = R2, nperm = 999), 
                           TRUE))
        if (class != "try-error") { 
          sign <- sort(fsel$order)
          select <- mem[, c(sign)] 
          list(MEM.select = select, NbVar = length(sign), pval = pval, R2.global = R2,
               R2.select = RsquareAdj(rda(response, select))$adj.r.squared, 
               AdjR2Cum = fsel$AdjR2Cum)
        } else return(pval)
      } else if (selection == "MIR") {
          mir <- MEM.moransel(as.vector(response)[, 1], attributes(mem)$listw, MEM.autocor = c, 
                              alpha = thresh)
          if (class(mir) == "list") {
            select <- mir$MEM.select
            list(MEM.select = select, NbVar = ncol(select), pval = pval, R2.global = R2, 
                 R2.select = RsquareAdj(rda(response, select))$adj.r.squared, 
                 AdjR2Cum = "Computed only for crit = forward")
          } else return(pval)
      } else if (selection == "global") {
          select <- mem
          list(MEM.select = select, NbVar = ncol(mem), pval = pval, R2.global = R2, 
               R2.select = R2, AdjR2Cum = "Computed only for crit = forward")
             }
      } else return(pval)
  }                                                  # End of the MEM.test() function
  # **********************************************************************************
  
  # Since the loop is entered only once if autocor = "positive" or "negative" 
  # but twice if autocor = "all" (once for the positive and once for the negative 
  # MEM models), we begin by setting up the loop.
  
  if (autocor != "all") { 
    k <- 1   # number of times the for loop will run
    if (autocor == "positive") { 
      cor <- c("positive", "negative") 
    }  else cor <- c("negative", "positive") 
  } else {
    k <- 2
    cor <- c("positive", "negative")
    # For the selection of the SWM on the basis of the sum of the pos. and neg. models:
    list.listW <- vector("list", 2)
    list.results <- vector("list", 2)
    list.listMEM <- vector("list", 2)
    list.listR2 <- vector("list", 2)
    
    # Internal function for the sum of the R2 of the pos. and neg. autocorrelated MEM variables:
    sum.R2 <- function(n) {
      # n is a vector of two values (R2 of the selected subset of MEM variables of a given SWM 
      # for its positive and negative MEM variables):
      nb.na <- length(which(is.na(n)) == TRUE)
      if (nb.na == 0) return(sum(n))
      else if (nb.na == 2) return(NA)
      else if (is.na(n[1]) == TRUE) return(n[2])
      else return(n[1])
    }
  }
  
  # A multitest p-value correction is needed for controling the type-I error rate. 
  # We define the total nb of tests:
  # ********************************
  # If one tests more than 1 SWM, 'candidates' only has one class, that is, 'list', and
  # the length of this list is the number of candidate W matrices. 
  # The same applies if ones generates a single SWM contained in a list (as it is done
  # with listw.candidates if only one SWM is generated), and if the argument 'candidates' of 
  # listw.select is this list. 
  # However, in the case of a single SWM, if 'candidates' is the listw object (and not a 
  # list containing it), then 'candidates' has two classes ('listw' and 'nb') and the length 
  # of 'candidates' does not correspond to the number of candidate SWMs anymore. 
  # This may occur if the SWM was not generated with listw.candidates:
  if (length(class(candidates)) == 1) {
    nbtest <- length(candidates)
    control <- TRUE
  } else {
    nbtest <- 1
    control <- FALSE
  }
  
  # 'struct' will be set to TRUE if at least one SWM is significant: 
  struct <- FALSE
  
  for (h in 1:k) {
    
    # For model comparison and selection:
    results <- as.data.frame(matrix(nrow = nbtest, ncol = 4))
    colnames(results) <- c("pval_sidak", "R2.global", "R2.select", "NbMEM")     
    # List of the MEM.select matrices (subsets of MEM variables for each significant SWM):
    listMEM <- vector("list", nbtest)
    # and corresponding R2:
    listR2 <- vector("list", nbtest)
    # List of 'MEM.test' results:
    listtest <- vector("list", nbtest)
    # List of global W matrices (all MEM variables):
    listW <- vector("list", nbtest)
    
    for (q in 1:nbtest) {
      if (control == TRUE) W <- scores.listw(candidates[[q]], MEM.autocor = cor[h], 
                                             store.listw = TRUE)
      else W <- scores.listw(candidates, MEM.autocor = cor[h], store.listw = TRUE)
      listW[[q]] <- W
      listtest[[q]] <- MEM.test(x, W, c = cor[h])
    }
    # Save the results in order to compare them and choose the best model:
    # ********************************************************************
    p <- c()   # To avoid signif. models for which no variable is selected by the forward sel.
    for (i in 1:nbtest) {
      if (is.list(listtest[[i]]) == TRUE) {
        results[i, 1] <- listtest[[i]]$pval
        results[i, 2] <- listtest[[i]]$R2.global
        results[i, 3] <- listtest[[i]]$R2.select
        results[i, 4] <- listtest[[i]]$NbVar
        
        listMEM[[i]] <- listtest[[i]]$MEM.select
        listR2[[i]] <- listtest[[i]]$AdjR2Cum
        p <- c(p, listtest[[i]]$pval)
      } else {
        results[i, 1] <- listtest[[i]]
        p <- c(p, 1)
      }
    }
    # If autocor == "all", then listW, results, listMEM, and listR2 must be saved in order not
    # to be erased by the second round of the for loop (when h = 2):
    if (autocor == "all") {
      list.listW[[h]] <- listW
      list.results[[h]] <- results
      list.listMEM[[h]] <- listMEM
      list.listR2[[h]] <- listR2
    }
    # Selection of the best SWM (and best eigenvector subset):
    # ********************************************************
    if (length(which(p <= alpha)) > 0) {
      struct <- TRUE
      if (autocor != "all") {
        if (crit == "forward") 
          best <- which.max(results[, 3])
        else if (crit == "global") 
          best <- which.max(results[, 2])
        else {   # crit = "MIR"
          rows_min <- which(results[, 4] == min(na.omit(results[, 4])))
          if (length(rows_min) == 1) best <- rows_min
          else best <- as.numeric(rownames(results)[rows_min][which.max(results[rows_min, 3])])
        }
      
        if (nbtest == 1) {
          if (control == TRUE) bestlistw = candidates[[1]]
          else bestlistw = candidates
        } else bestlistw = candidates[[best]]
      
        L <- list(MEM.all = listW[[best]], MEM.select = listMEM[[best]], listw = bestlistw,
                  MEM.AdjR2Cum = listR2[[best]], name = names(candidates)[best], 
                  pval = results[best, 1], R2.global = results[best, 2], 
                  R2.select = results[best, 3], NbVar = results[best, 4], bestw_index = best)
        L.all <- list(MEM.all = listW, pval = results[, 1], R2.global = results[, 2], 
                      R2.select = results[, 3])
      }
    }
  }    # End of the loop 'for (h in 1:k)'
  
  # if autocor = "all" and if a signif. pattern was detected at least once for the positive and 
  # the negative MEM variables, then the best model is computed from the sum of both. However, 
  # it is not because both a positive and a negative MEM model have been detected at least once
  # that the best SWM will contain both:
  if (autocor == "all" & struct == TRUE) {
    if (crit == "forward") {
      bind <- cbind(list.results[[1]][, 3], list.results[[2]][, 3])
      sum <- apply(bind, 1, sum.R2)
      best <- which.max(sum)
    } else if (crit == "global") {
      bind <- cbind(list.results[[1]][, 2], list.results[[2]][, 2])
      sum <- apply(bind, 1, sum.R2)
      best <- which.max(sum)
    } else { # crit = "MIR"
      bind <- cbind(list.results[[1]][, 4], list.results[[2]][, 4])                        
      sum <- apply(bind, 1, sum.R2)                       
      rows_min <- which(sum == min(na.omit(sum)))
      if (length(rows_min) == 1) 
        best <- rows_min
      else {
        bind.R2 <- cbind(list.results[[1]][, 3], list.results[[2]][, 3])
        sumR2 <- apply(bind.R2, 1, sum.R2)
        best <- c(1:length(sum))[rows_min][which.max(sumR2[rows_min])]
      }
    }
    
    if (nbtest == 1) {
      if (control == TRUE) bestlistw = candidates[[1]]
      else bestlistw = candidates
    } else bestlistw = candidates[[best]]

    # To recover the global models containing all positive and negative MEM variables together:
    listW2 <- vector("list", length(listW))    
    for (i in 1: nbtest) {
      nb.pos <- ncol(list.listW[[1]][[i]])
      nb.neg <- ncol(list.listW[[2]][[i]])
      colnames(list.listW[[2]][[i]]) <- paste("MEM", seq(nb.pos + 1, nb.pos + nb.neg), sep = "")
      listW2[[i]] <- cbind(list.listW[[1]][[i]], list.listW[[2]][[i]])
      attributes(listW2[[i]])$values <- c(attributes(list.listW[[1]][[i]])$values, 
                                          attributes(list.listW[[2]][[i]])$values)
    }
    # To join the selected MEM variables within the pos. and neg. portions of the best SWM:
    if (crit == "global") {
      listMEM <- listW2
      listR2 <- "Computed only for crit = forward"
    } else {
      # We check whether the best SWM contains both positive and negative MEM variables.
      # If it does, then we correct the ID (number) of the negative MEM variables (they were
      # computed separately, so that the first negative MEM variable is wrongly called "MEM1"):
      if (length(is.na(c(list.results[[1]][best, 3], list.results[[2]][best, 3]))) == 0) {
        nb.pos <- ncol(list.listW[[1]][[best]])
        col <- colnames(list.listMEM[[2]][[best]])
        id.neg <- as.numeric(substr(col, 4, nchar(col[length(col)]))) + nb.pos
        colnames(list.listMEM[[2]][[best]]) <- paste("MEM", id.neg, sep = "")
        listMEM <- cbind(list.listMEM[[1]][[best]], list.listMEM[[2]][[best]])
        if (crit == "forward") 
          listR2 <- c(list.listR2[[1]][[best]], list.listR2[[2]][[best]])
        else listR2 <- "Computed only for crit = forward"
      } else if (is.na(list.results[[1]][best, 3])) {
        nb.pos <- ncol(list.listW[[1]][[best]])
        col <- colnames(list.listMEM[[2]][[best]])
        id.neg <- as.numeric(substr(col, 4, nchar(col[length(col)]))) + nb.pos
        colnames(list.listMEM[[2]][[best]]) <- paste("MEM", id.neg, sep = "")
        listMEM <- list.listMEM[[2]][[best]]
        if (crit == "forward") 
          listR2 <- list.listR2[[2]][[best]]
        else listR2 <- "Computed only for crit = forward"
      } else {
        listMEM <- list.listMEM[[1]][[best]]
        if (crit == "forward") 
          listR2 <- list.listR2[[1]][[best]]
        else listR2 <- "Computed only for crit = forward"
      }
    }
     
    L <- list(MEM.all = listW2[[best]], MEM.select = listMEM, listw = bestlistw,
              MEM.AdjR2Cum = listR2, name = names(candidates)[best], 
              pval = c(MEM_positive = list.results[[1]][best, 1], 
                       MEM_negative = list.results[[2]][best, 1]), 
              R2.global = c(MEM_positive = list.results[[1]][best, 2], 
                            MEM_negative = list.results[[2]][best, 2]), 
              R2.select = c(MEM_positive = list.results[[1]][best, 3], 
                            MEM_negative = list.results[[2]][best, 3]), 
              NbVar = c(MEM_positive = list.results[[1]][best, 4], 
                        MEM_negative = list.results[[2]][best, 4]), bestw_index = best)
    L.all <- list(MEM.all = listW2, pval = cbind(MEM_positive = list.results[[1]][, 1], 
                                                 MEM_negative = list.results[[2]][, 1]), 
                  R2.global = cbind(MEM_positive = list.results[[1]][, 2], 
                                    MEM_negative = list.results[[2]][, 2]),
                  R2.select = cbind(MEM_positive = list.results[[1]][, 3], 
                                    MEM_negative = list.results[[2]][, 3]))
  }   # End of the loop 'if (autocor == "all" & struct == TRUE)'
  
  # Output of the function:
  # ***********************
  if (autocor == "all" & struct == TRUE) {
    if (summary == TRUE) {
      if (length(which(L$pval <= alpha)) == 2) {
        if (crit != "global") val <- c(round(L$R2.select[1], 3), round(L$R2.select[2], 3))
        else val <- c(round(L$R2.global[1], 3), round(L$R2.global[2], 3))
      cat("\n", 
          "***********************************************************************************",
          "\n",
          "A best model containing both positive (", sidak, " p-value = ", round(L$pval[1], 5),
          " , adjusted R2 = ", val[1], ")", "\n", "and negative MEM variables (", sidak, 
          " p-value = ", round(L$pval[2], 5), ", adjusted R2 = ", val[2], 
          ") was selected.", "\n", 
          "***********************************************************************************",
          "\n", sep = "")
      } else {
        nb <- which(L$pval <= alpha)
        sign <- cor[nb]
        if (crit != "global") val <- round(L$R2.select[nb], 3)
        else val <- round(L$R2.global[nb], 3)
        cat("\n", 
            "*******************************************************************************",
            "****", "\n",
            "A best model containing only ", sign, " MEM variables was selected (", sidak, 
            " p-value", "\n", "= ", round(L$pval[nb], 5), ", adjusted R2 = ", val, ").", "\n", 
            "*******************************************************************************",
            "****", "\n", sep = "")
      }
    }
    return(list(all = L.all, best = L))
  } else if (struct) {
      if (summary == TRUE) {
        if (crit != "global") val <- round(L$R2.select, 3)
        else val <- round(L$R2.global, 3)
        cat("\n", 
            "*******************************************************************************",
            "****", "\n", 
            "A best ", cor[1], " MEM model was selected (", sidak, " p-value = ", 
            round(L$pval, 5), ", adjusted R2 = ", val, ").", "\n", 
            "*******************************************************************************", 
            "****", "\n", sep = "")
      }
      return(list(all = L.all, best = L))
  } else {
      if (summary == TRUE)
        cat("\n", "**********************************************************", 
            "\n", "No significant spatial structure was detected in the data.", "\n", sep = "")
  }
}

#' @rdname listw.select
"mem.select" <- function(x, 
                         candidates, 
                         autocor = c("positive", "negative", "all"), 
                         crit = c("forward", "MIR", "global"),
                         alpha = 0.05, 
                         summary = TRUE) {
  crit <- match.arg(crit)
  autocor <- match.arg(autocor)
  res <- listw.select(x = x, candidates = candidates, autocor = autocor, crit = crit,
                      alpha = alpha, correction = TRUE, summary = summary)
  res2 <- list(MEM.all = res$best$MEM.all, MEM.select = res$best$MEM.select, 
               listw = res$best$listw, MEM.AdjR2Cum = res$best$MEM.AdjR2Cum, 
               pval = res$best$pval, R2.global = res$best$R2.global, 
               R2.select = res$best$R2.select, NbVar = res$best$NbVar)
  return(res2)
}
  