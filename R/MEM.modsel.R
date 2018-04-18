#' Function to optimize the selection of a spatial weighting matrix and select the best 
#' subset of eigenvectors
#' 
#' \code{MEM.modsel} computes spatial eigenvectors for various definitions of spatial weighting
#' matrices (W matrices) and optimizes the selection of the W matrix by maximizing the
#' value of the adjusted R-squared (R2) while controling the type I error rate (see
#' Bauman et al. 2018a). The function also selects the best subset of eigenvectors to be
#' used as spatial predictors within the best W matrix by performing a forward selection with
#' double stopping criterion (Blanchet et al. 2008). The best W matrix returned is selected
#' on the basis of the R2 of either the global models or the subsets of eigenvectors selected
#' by the forward selections.
#' \code{MEM.modsel} combines calls to the functions \code{scores.listw} and \code{forward.sel}.
#' The list of candidate W matrices can easily be generated using the user-friendly 
#' \code{listw.candidates} function. The significance of each W matrix is tested by
#' 9999 permutations with the function \code{anova.cca} (package \code{vegan}).
#' 
#' @details While the selection of the W matrix is the most critical step of the spatial
#' eigenvector-based methods (Dray et al. 2006), Bauman et al. (2018a) showed that 
#' optimizing the choice of the W matrix led to inflated type I error rates if an
#' explicit control of the number of W matrices tested was not applied. The function
#' \code{MEM.modsel} therefore applies a Sidak correction (Sidak 1967) for multiple tests to
#' the p-value of the global test of each W matrix (i.e., the model integrating the 
#' whole set of spatial predictors). The Sidak correction is computed as: 
#'\eqn{P_corrected = 1 - (1 - P)^n}, where \eqn{n} is the number of tests performed, \eqn{P}
#' is the observed p-value, and \eqn{P_corrected} is the new p-value after the correction.
#' The p-value is first computed by 9999 permutations with the function \code{anova.cca} 
#' (package \code{vegan}), and is then corrected according to the total number of W matrices 
#' tested (if \code{correction = TRUE}). Although the function can be run without this
#' correction (\code{correction = FALSE}), using the default value (\code{correction = TRUE})
#' is strongly recommended to avoid inflated type I error rates (Bauman et al. 2018a).
#' 
#' As a consequence of the necessity of the p-value correction, the significance threshold
#' decreases as the number of W matrices tested increases, hence leading to a trade-off
#' between the gain of accuracy and the statistical power loss. See \code{Details} of the 
#' function \code{\link{listw.candidates}} for recommendations about the selection of a 
#' suitable set of candidate W matrices. 
#' 
#' The optimization criterion of the W matrix performed by \code{MEM.modsel} is either based 
#' on all the positively or negatively autocorrelated eigenvectors (\code{crit = "global"}), 
#' or is based on an optimized subset of spatial eigenvectors (also referred to as spatial 
#' predictors or MEM variables) (\code{crit = "forward"} or \code{"MIR"}).
#' If one only wants to optimize the selection of the W matrix, without the intervention of the
#' selection of a subset of predictors within each W matrix (\code{crit = "global"}), then the
#' best W matrix returned will be selected among the significant W matrices as the one for which
#' the complete set of (positively or negatively autocorrelated) MEM variables yields the
#' highest adjusted R2. This may be interesting when actually a global set of MEM variables,
#' like in Moran spectral randomizations (Wagner and Dray 2015) or when using smoothed MEM
#' (Munoz 2009). If the MEM variables will be further used in a model including 
#' actual predictors (e.g. environmental), then a subset of spatial eigenvectors will need to
#' be selected before proceeding to further analyses in order to avoid model overfitting 
#' and a loss of statistical power to detect the contribution of the environment to the 
#' variability of the response data (Griffith 2003, Dray et al. 2006, Blanchet et al. 2008,
#' Peres-Neto and Legendre 2010, Diniz-Filho et al. 2012). In such cases, the optimization
#' of the W matrix should be based on the best subset of spatial predictors within each
#' candidate W matrix. Although several eigenvector selection approaches have been proposed to 
#' select this best subset of eigenvectors, Bauman et al. (2018b) showed that two main 
#' procedures should be preferred, depending on the underlying objective.
#' 
#' The most powerful and accurate selection method, in terms of R2 estimation, is the 
#' forward selection with double stopping criterion of Blanchet et al. (2008). This method
#' (\code{crit = "forward"}, default option) should be preferred when the objective is to
#' describe as accurately as possible the spatial patterns of the response data (\code{x}).
#' If the objective is to optimize the detection of the spatial patterns in the residuals of
#' a model of the response variable(s) against a set of environmental predictors, for instance,
#' then \code{x} should be the model residuals, and \code{crit = "forward"} or \code{"MIR"}.
#' The latter option would be more suitable if one wants to selects the smallest subset of 
#' spatial predictors necessary to capture all the spatial autocorrelation in \code{x} (see
#' function \code{\link{MEM.moransel}}. This criterion has the advantage to maintain the 
#' standard errors of the actual predictor coefficients as low as possible. 
#' If the objective is to describe the spatial patterns of \code{x} as accurately as possible, 
#' then the former option is more appropriate (but is likely to select a slightly higher 
#' number of spatial predictors, see Bauman et al. 2018b).
#' If \code{crit = "forward"}, \code{MEM.modsel} performs the forward selection on the 
#' significant W matrices and selects among these the W matrix for which the forward-selected
#' subset of spatial eigenvectors yields the highest adjusted R2. If \code{crit = "MIR"},
#' \code{MEM.modsel} performs the MIR selection on all the significant candidate W matrices,
#' and selects the best W matrix as the one with the smallest number of MIR-selected spatial 
#' eigenvectors. If two or more W matrices present the same smallest number of predictors, 
#' then the choice is made among them on the basis of the adjusted R2, as for a same number of 
#' spatial eigenvectors, a higher R2 indicate a more accurate description of the spatial 
#' structure in the residuals (Bauman et al. 2018a).
#'  
#' @param x Vector, matrix, or dataframe of response data
#' @param candidates A list of one or more spatial weighting matrices of the class 
#' \code{listw}; \code{candidates} can be a list containing one or several spatial 
#' weighting matrices, but it can also be a single listw object.
#' @param autocor Sign of the spatial eigenvectors to consider; "positive", "negative", 
#' or "all", for positively, negatively autocorrelated eigenvectors, or both, respectively; 
#' default is "positive"
#' @param crit Criterion to select the best W matrix among the set of significant candidate 
#' matrices. Either \code{forward} (default option), \code{"MIR"}, or \code{"global"} (see
#' \code{Details}).
#' @param alpha_thresh Uncorrected predefined significance value; default is 0.05
#' @param correction wheter the p-value of the global test performed on each W matrix should
#' be corrected for multiple tests (TRUE) or not (FALSE); default is TRUE
#' @param summary Logical; Whether a message summarizing the results should be returned
#' or not; The R2 value returned in the message depends on
#' the argument \code{crit} (R2 of the subset of forward selected predictors or R2 of the
#' global model); Default is \code{TRUE}
#' 
#' @return The function returns two lists containing several features of the W matrix selected
#' by the optimization procedure, and general features of all the W matrices generated and 
#' tested: 
#' The first list, \code{all}, contains: 
#' 
#' \describe{
#' \item{MEM.all}{A list of all the total set of MEM variables (as generated by 
#' \code{scores.listw}) for the candidate W matrices (\code{listw} objects) provided to the 
#' function, in the same order as the W matrices in \code{candidates}.}
#' \item{pval}{The p-values (corrected by the Sidak correction if \code{correction = TRUE})
#' of the global models for each W matrix.}
#' \item{R2.global}{The adjusted R2 of the complete set of MEM variables of each W 
#' matrix. R2.global is available only for the significant W matrices.}
#' }
#' 
#' The second list, \code{best}, contains the features of the best model: 
#' 
#' \describe{
#' \item{MEM.all}{Dataframe of the complete set of eigenvectors of the best W matrix.}
#' \item{MEM.select}{Dataframe of the subset of eigenvectors selected by the criterion
#' \code{crit}. This set of spatial predictors is the one to be used in further analyses 
#' (e.g. variation partitioning, ordinary least squares, GLM) (see \code{Details} section). 
#' If \code{crit = "global"}, then \code{MEM.select} contains the entire set of MEM variables 
#' generated for the corresponding W matrix. In this case, \code{MEM.select} and \code{MEM.all}
#' are the same.} 
#' \item{listw}{Element of \code{candidates} corresponding to the best W matrix. It is an 
#' object of class \code{listw}.}
#' \item{MEM.AdjR2Cum}{Only for \code{crit = "forward"}; Vector of the cumulative adjusted 
#' R2 of the eigenvectors selected by the forward selection (\code{MEM.select}).}
#' \item{name}{Name of the best W matrix ("Connectivity matrix_Weighting matrix") if the 
#' list of candidate W matrices was generated by the function \code{\link{listw.candidates}}. 
#' Otherwise, corresponds to the index of the best W matrix in the list \code{candidates} 
#' provided to \code{MEM.modsel}.} 
#' \item{pval}{Global p-value of the best W matrix (after multiple-test correction, if 
#' the \code{correction} argument = TRUE).}
#' \item{R2.global}{The adjusted R2 of the complete set of MEM variables of the best W 
#' matrix.}
#' \item{R2.select}{Adjusted R2 of \code{MEM.select}. If \code{crit = "global"}, then 
#' \code{R2.select} is equal to \code{R2.global}.}
#' \item{NbVar}{Number of spatial predictors selected by the forward selection.}
#' \item{bestw_index}{Index of the best W matrix in the list \code{candidates} provided 
#' to \code{MEM.modsel}.}
#' }
#' 
#' If \code{autocor = "all"} and a significant spatial model was selected separately for the
#' MEM variables associated to positive and negative eigenvalues, then both lists
#' (\code{all} and \code{best}) are generated twice, for both sets of MEM variables,
#' and are returned in two different lists: 
#' 
#' \describe{
#' \item{MEM.pos}{Lists \code{all} and \code{best} for the MEM variables corresponding to
#' a positive autocorrelation.}
#' \item{MEM.neg}{Lists \code{all} and \code{best} for the MEM variables corresponding to
#' a negative autocorrelation.}
#' }
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
#' with irregular sampling. Ecological Modelling, 220, 2683–-2689.
#' 
#' Peres-Neto P. and Legendre P. (2010) Estimating and controlling for spatial structure 
#' in the study of ecological communities. Global Ecology and Biogeography, 19, 174--184
#' 
#' Sidak Z. (1967) Rectangular confidence regions for the means of 
#' multivariate normal distributions. Journal of the American Statistical Association, 62(318),
#' 626--633
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
#' ### we want to test and compare (with the MEM.modsel function). We test a Gabriel's graph, 
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
#' ### Optimization of the selection of the W matrix among the candidates generated above, 
#' ### using the corrected significance threshold calculated above for the global tests:
#' W_sel <- MEM.modsel(Y_sampled, candidates, autocor = "positive", crit = "forward", 
#'                     correction = TRUE)
#' ### Some characteristics of the best spatial model:
#' # Best W matrix:
#' W_sel$best$name
#' # Selected subset of spatial predictor within the best W matrix:
#' W_sel$best$MEM.select
#' W_sel$best$NbVar
#' # Corrected p-value of the global test of the best W matrix: 
#' W_sel$best$pval
#' # Adjusted R2 of the subset of spatial predictors selected within the chosen W matrix:
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
#' @importFrom vegan rda anova.cca RsquareAdj
#' @export MEM.modsel

"MEM.modsel" <- function(x, 
                         candidates, 
                         autocor = c("positive", "negative", "all"), 
                         crit = c("forward", "MIR", "global"),
                         alpha_thresh = 0.05, 
                         correction = TRUE,
                         summary = TRUE) {
  
  x <- as.data.frame(x)
  if (any(is.na(x))) stop ("NA entries in x")
  
  if (correction == TRUE) sidak <- "corrected" else sidak <- "uncorrected"
  autocor <- match.arg(autocor)
  
  crit <- match.arg(crit)

  # **********************************************************************************  
  MEM.test <- function(response = x, 
                       mem, 
                       selection = crit,
                       c = autocor, 
                       d = nbtest, 
                       alpha = alpha_thresh,
                       corr = correction) {
    # MEM.test is an internal function that tests the significance of a W matrix while taking
    # into consideration the total number of W matrices tested in MEM.modsel (if corr = TRUE).
    # This total number of tests is used to apply a correction to the p-value in order not to
    # inflate the type I error rate (if corr = TRUE). If the tested W matrix is significant 
    # at the corrected threshold value of significance, then either a model selection is 
    # performed using Blanchet et al.'s forward selection (2008) (crit = "forward") or the
    # MIR selection (crit = "MIR"), or only the adjusted R2 of the global model is computed
    # (crit = "global").
    pval <- anova.cca(rda(response, mem), permutations = 9999)$Pr[1]
    if (correction == TRUE) pval <- 1-(1-pval)^d    # Sidak correction 
    if (c == "all") pval <- 1-(1-pval)^2            # Sidak correction (autocor= "all") 
    if (pval <= alpha) {  
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
          mir <- MEM.moransel(response, attributes(mem)$listw, MEM.autocor = c, alpha = alpha)
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
  }
  
  # A multitest p-value correction is needed for controling the type-I error rate. 
  # We define the total nb of tests:
  # ********************************
  # If one tests more than 1 W matrix, 'candidates' only has one class, that is, 'list', and
  # the length of this list is the number of candidate W matrices. 
  # The same applies if ones generates a single W matrix contained in a list (as it is done
  # with listw.candidates if only one W matrix is generated), and if the argument 'candidates' 
  # is this list. However, in the case of a single W matrix, if 'candidates' is the listw 
  # object (and not a list containing it), then 'candidates' has two classes ('listw' and 'nb')
  # and the length of 'candidates' does not correspond to the number of candidate W matrices
  # anymore. This may occur if the W matrix was not generated with listw.candidates:
  if (length(class(candidates)) == 1) {
    nbtest <- length(candidates)
    control <- TRUE
  } else {
    nbtest <- 1
    control <- FALSE
  }
  
  lenlist <- c()   # Will help with the result output
  
  for (h in 1:k) {
    
    # For model comparison and selection:
    results <- as.data.frame(matrix(nrow = nbtest, ncol = 4))
    colnames(results) <- c("pval_sidak", "R2.global", "R2.select", "NbMEM")     
    # List of the MEM.select matrices (subsets of MEM variables for each signif. W matrix):
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
      listtest[[q]] <- MEM.test(x, W)
    }
    # Save the results in order to compare them and choose the best model:
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
    # Selection of the best W matrix (and best eigenvector subset within it):
    if (length(which(p <= alpha_thresh)) > 0) {
      
      if (crit == "forward") best <- which.max(results[, 3])
      else if (crit == "global") best <- which.max(results[, 2])
      else {   # crit = "MIR"
        min_Nb_MEM <- which(results[, 4] == which.min(results[, 4]))
        if (length(min_Nb_MEM) == 1) best <- which.min(results[, 4])
        else best <- which.max(results[, 3][min_Nb_MEM])
      }
      
      if (nbtest == 1) {
        if (control == TRUE) bestlistw = candidates[[1]]
        else bestlistw = candidates
      } else bestlistw = candidates[[best]]
      if (autocor != "all") {
        lenlist <- c(lenlist, cor[h])
        L <- list(MEM.all = listW[[best]], MEM.select = listMEM[[best]], listw = bestlistw,
                  MEM.AdjR2Cum = listR2[[best]], name = names(candidates)[best], 
                  pval = results[best, 1], R2.global = results[best, 2], 
                  R2.select = results[best, 3], NbVar = results[best, 4], bestw_index = best)
        L.all <- list(MEM.all = listW, pval = results[, 1], R2.global = results[, 2], 
                      R2.select = results[, 3])
      } else {
        lenlist <- c(lenlist, cor[h])
        if (h == 1) {
          L1 <- list(MEM.all = listW[[best]], MEM.select = listMEM[[best]], 
                     listw = bestlistw, MEM.AdjR2Cum = listR2[[best]], 
                     name = names(candidates)[best], pval = results[best, 1], 
                     R2.global = results[best, 2], R2.select = results[best, 3], 
                     NbVar = results[best, 4], bestw_index = best)
          L1.all <- list(MEM.all = listW, pval = results[, 1], R2.global = results[, 2],
                         R2.select = results[, 3])
        } else {
          L2 <- list(MEM.all = listW[[best]], MEM.select = listMEM[[best]], 
                     listw = bestlistw, MEM.AdjR2Cum = listR2[[best]], 
                     name = names(candidates)[best], pval = results[best, 1], 
                     R2.global = results[best, 2], R2.select = results[best, 3], 
                     NbVar = results[best, 4], bestw_index = best)
          L2.all <- list(MEM.all = listW, pval = results[, 1], R2.global = results[, 2],
                         R2.select = results[, 3])
        }
      }
    }
  }    # End of the for loop
  
  # Output of the MEM.modsel function:
  if (length(lenlist) == 2) {
    if (summary == TRUE) {
      if (crit != "global")
        val <- c(round(L1$R2.select, 3), round(L2$R2.select, 3))
      else 
        val <- c(round(L1$R2.global, 3), round(L2$R2.global, 3))
      cat("\n", "*****************************************************", "\n",
          "*****************************************************", "\n",
          "A best positive (", sidak, " p-value = ", round(L1$pval, 5), ", adjusted R2 =", 
          "\n", val[1], ") and best negative (", sidak, " p-value = ", round(L2$pval, 5), ",", 
          "\n", "adjusted R2 = ", val[2], ")", "MEM models were selected.", "\n", 
          "*****************************************************", "\n",
          "*****************************************************", "\n",
          "The output of the function is a list of two lists (MEM.pos and MEM.neg).", "\n",
          sep = "")
    }
    return(list(MEM.pos = list(all = L1.all, best = L1), 
                MEM.neg = list(all = L2.all, best = L2)))
  } else 
    if (length(lenlist) == 1) { 
      if (autocor == "all") {
        if (lenlist == "positive") {
          L <- L1
          L.all <- L1.all
        } else {
          L <- L2
          L.all <- L2.all
        }
      }
      if (summary == TRUE) {
        if (crit != "global")
          val <- round(L$R2.select, 3)
        else 
          val <- round(L$R2.global, 3)
        cat("\n", "*****************************************************", 
            "\n", "*****************************************************", "\n",
            "A best ", lenlist, " MEM model was selected (", sidak, " p-value = ", 
            round(L$pval, 5), ", adjusted R2 = ", val, ").", "\n", 
            "*****************************************************", "\n",
            "*****************************************************", "\n", sep = "")
      }
      return(list(all = L.all, 
                  best = list(MEM.all = L$MEM.all, MEM.select = L$MEM.select, 
                              listw = L$listw, MEM.AdjR2Cum = L$MEM.AdjR2Cum, 
                              name = L$name, pval = L$pval, R2.global = L$R2.global, 
                              R2.select = L$R2.select, NbVar = L$NbVar, 
                              bestw_index = L$bestw_index)))
    } else 
      if (summary == TRUE)
        cat("\n", "*****************************************************", 
            "\n", "*****************************************************", "\n",
            "No significant spatial structure was detected in the data.", "\n",
            "\n", sep = "")
}
