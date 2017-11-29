#' Function to optimise the selection of a spatial weighting matrix and select the best 
#' subset of eigenvectors within it
#' 
#' MEM.modsel computes spatial eigenvectors for various definitions of spatial weighting
#' matrices (W matrices) and optimises the selection of the W matrix by maximising the
#' value of the adjusted R-square while controling the type I error rate (see
#' Bauman et al. 2018). The function also selects the best subset of eigenvectors to be
#' used as spatial predictors within the best W matrix by performing a forward selection.
#' It combines calls to the functions \code{scores.listw} and \code{forward.sel}.
#' The list of W matrix candidates can easily be generated using the user-friendly 
#' \code{listw.candidates} function. The significance of each W matrix is tested by
#' 9999 permutations by the means of the function \code{anova.cca} (package \code{vegan}).
#' 
#' @details While the selection of the W matrix is the most critical step of the spatial
#' eigenvector-based methods (Dray et al. 2006), Bauman et al. (2018) showed that 
#' optimising the choice of the W matrix led to inflated type I error rates if an
#' explicit control of the number of W matrices tested was not applied. The function
#' MEM.modsel therefore applies a Sidak correction (Sidak 1967) for multiple tests to
#' the p-value of the global test of each W matrix (i.e., the model integrating the 
#' whole set of spatial predictors). The Sidak correction is computed as: 
#'\eqn{P_corrected = 1 - (1 - P)^n}, where \eqn{n} is the number of tests performed, \eqn{P}.
#' is the observed p-value, and \eqn{P_corrected} is the new p-value after the correction.
#' The p-value is first computed by 9999 permutations with the function \code{anova.cca} 
#' (package \code{vegan}), and is then corrected according to the total number of W matrices 
#' tested (if \code{correction = TRUE}). Although the function can be run without this
#' correction (\code{correction = FALSE}), using the default value (\code{correction = TRUE})
#' is strongly recommended to avoid highly inflated type I error rates (Bauman et al. 2018).
#' As a consequence of the necessity of the p-value correction, the significance threshold
#' decreases as the number of W matrices tested increases, hence leading to a trade-off
#' between the gain of accuracy and the statistical power loss. See 'Details' of the 
#' function \code{\link{listw.candidates}} for recommendations about selecting a suitable set 
#' of W matrix candidates. 
#' Once the W matrix was optimised, a subset of eigenvectors should be selected before 
#' proceding to further analyses in order both to avoid model overfitting and a loss of 
#' statistical power to detect the contribution of the environment to the variability of 
#' the response data (Griffith 2003, Dray et al. 2006, Blanchet et al. 2008, 
#' Peres-Neto and Legendre 2010, Diniz-Filho et al. 2012). Although several eigenvector
#' selection approaches have been proposed, Bauman et al. (2017) showed that the most
#' powerful and accurate method, in terms of adjusted R-square estimation, was the 
#' forward selection with double stopping criterion of Blanchet et al. (2008). 
#' The function \code{MEM.modsel} performs the forward selection on the significant W 
#' matrices and selects among these the best W matrix as the one for which the 
#' eigenvector subset resulting from the forward selection displays the highest adjusted R².
#'  
#' @param x Vector, matrix, or dataframe of response data
#' @param candidates A list of one or more spatial weighting matrices of class listw
#' @param autocor Sign of the spatial eigenvectors to consider; "positive", "negative", 
#' or "all", for positively, negatively autocorrelated eigenvectors, or both, respectively; 
#' default is "positive"
#' @param alpha_thresh Uncorrected predefined significance value; default is 0.05
#' @param correction wheter the p-value of the global test performed on each W matrix should
#' be corrected for multiple tests (TRUE) or not (FALSE); default is TRUE
#' 
#' @return The function returns a list containing several features of the W matrix selected
#' by the optimisation procedure: \describe{
#' \item{MEM.all }{Dataframe of the complete set of eigenvectors of the best W matrix.} 
#' \item{MEM.select }{Dataframe of the eigenvectors selected by the forward selection. These
#' are the spatial predictors that should be used for further analyses. In the paradigm of
#' 'spatial legacy' (Peres-Neto and Legendre 2010), these spatial predictors are used to
#' model the spatial patterns of the response data as accurately as possible. For instance,
#' using them together with a set of environmental variables in a variation partitioning 
#' (see function \code{\link{varpart}} from the \code{vegan} package) allows assessing the
#' pure and joint contribution of the environmental and spatial predictors. In the 
#' 'spatial nuisance' paradigm (Peres-Neto and Legendre 2010), the spatial predictors are
#' used to remove the spatial signal from the response data in order to respect the condition
#' of independency of the residuals in classical OLS or GLM models (spatial filtering through
#' spatial eigenvector mapping for instance, Diniz-Filho et al. 2005). See Bauman et al.
#' (2017) for a discussion on the choice of the spatial eigenvector selection method 
#' according to the univariate or multivariate nature of the response data, and to the
#' spatial paradigm considered.} 
#' \item{listw}{Element of \code{candidates} (and object of class \code{listw}) corresponding 
#' to the best W matrix selected.}
#' \item{MEM.AdjR2Cum }{Vector of the cumulative adjusted R² of the eigenvectors of 
#' \code{MEM.select }.} 
#' \item{name }{Name of the best W matrix ("Connectivity matrix_Weighting matrix") if the 
#' list of W matrix candidates was generated by the function \code{\link{listw.candidates}}. 
#' Otherwise, corresponds to the index of the best W matrix in the list \code{candidates} 
#' provided to \code{MEM.modsel}.} 
#' \item{pval }{Global p-value of the best W matrix (after multiple-test correction, if 
#' the \code{correction} argument = TRUE).} 
#' \item{R2adj }{Adjusted R-square of \code{MEM.select}.} 
#' \item{NbVar }{Number of spatial predictors selected by the forward selection.} 
#' \item{bestw_index }{Index of the best W matrix in the list \code{candidates} provided 
#' to \code{MEM.modsel}.}}
#' 
#' @author Bauman David \email{dbauman@@ulb.ac.be} or \email{davbauman@@gmail.com}
#' 
#' @seealso \code{\link{listw.candidates}}, \code{\link{scores.listw}}, 
#' \code{\link[vegan]{varpart}}
#' 
#' @references Sidak Z. (1967) Rectangular confidence regions for the means of 
#' multivariate normal distributions. Journal of the American Statistical Association, 62(318),
#' 626--633
#' 
#' Griffith D. (2003) Spatial autocorrelation and spatial filtering: gaining understanding 
#' through theory and scientific visualization. Springer, Berlin
#' 
#' Dray S., Legendre P. and Peres-Neto P. R. (2006) Spatial modeling: a comprehensive 
#' framework for principal coordinate analysis of neighbor matrices (PCNM). Ecological 
#' Modelling, 196, 483--493
#' 
#' Blanchet G., Legendre P. and Borcard D. (2008) Forward selection of explanatory variables.
#' Ecology, 89(9), 2623--2632
#' 
#' Peres-Neto P. and Legendre P. (2010) Estimating and controlling for spatial structure 
#' in the study of ecological communities. Global Ecology and Biogeography, 19, 174--184
#' 
#' Diniz-Filho J.A.F., Bini L.M., Rangel T.F., Morales-Castilla I. et al. (2012) On the 
#' selection of phylogenetic eigenvectors for ecological analyses. Ecography, 35, 239--249
#' 
#' Bauman D., Drouet T., Dray S. and Vleminckx J. (2017) Disentangling good from bad 
#' practices in the selection of spatial or phylogenetic eigenvectors. Ecography 
#' *******************************
#' 
#' Bauman D., Fortin M-J, Suez M., Drouet T. and Dray S. (2018) *************************
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
#' ### the defaut value of DBthresh; see help of function listw.candidates). 
#' ### These connectivity matrices are then either not weighted (binary weighting), or 
#' ### weighted by the linearly decreasing function (see help of the function listw.candidates):
#' candidates <- listw.candidates(coord = xy, del = FALSE, rel = FALSE, fconcdown = FALSE, 
#'                                fconcup = FALSE, PCNM = FALSE)
#' ### Number of W matrix candidates generated:
#' nbw <- length(candidates)
#' ### Significance threshold value after p-value correction (Sidak correction):
#' 1 - (1 - 0.05)^(1/nbw)
#' ### Optimisation of the selection of the W matrix among the candidates generated above, 
#' ### using the corrected significance threshold calculated above for the global tests:
#' W_sel <- MEM.modsel(Y_sampled, candidates, autocor = "positive", correction = TRUE)
#' W_sel$name
#' W_sel$MEM.select
#' W_sel$pval
#' W_sel$R2adj
#' }
#'  
#' @importFrom vegan rda anova.cca RsquareAdj
#' @export MEM.modsel

"MEM.modsel" <- function(x, 
                         candidates, 
                         autocor = c("positive", "negative", "all"), 
                         alpha_thresh = 0.05, 
                         correction = TRUE) {
  
  x <- as.data.frame(x)
  if (any(is.na(x))) stop ("NA entries in x")
  
  if (correction == TRUE) sidak <- "corrected" else sidak <- "uncorrected"
  autocor <- match.arg(autocor)

  # **********************************************************************************  
  MEM.test <- function(a = x, 
                       b, 
                       c = autocor, 
                       d = nbtest, 
                       alpha = alpha_thresh,
                       corr = correction) {
    # MEM.test is an internal function that tests the significance of a W matrix while taking
    # into consideration the total number of W matrices tested in MEM.modsel. This total
    # number of tests is used to apply a correction to the p-value in order not to
    # inflate the type I error rate. If the tested W matrix is significant at the corrected
    # threshold value of significance, a model selection is performed using Blanchet et 
    # al.'s (2008) forward selection with double stopping criterion.
    pval <- anova.cca(rda(a, b), permutations = 9999)$Pr[1]
    if (correction == TRUE) pval <- 1-(1-pval)^d    # Sidak correction 
    if (c == "all") pval <- 1-(1-pval)^2            # Sidak correction (autocor= "all") 
    if (pval <= alpha) {  
      R2adj <- RsquareAdj(rda(a, b))$adj.r.squared
      class <- class(try(fsel <- forward.sel(a, b, adjR2thresh = R2adj, nperm = 999), TRUE))
      if (class != "try-error") { 
        sign <- sort(fsel$order)
        MEM.select <- b[, c(sign)] 
        list(MEM.select = MEM.select, NbVar = length(sign), pval = pval, 
             R2adj = RsquareAdj(rda(a, MEM.select))$adj.r.squared, AdjR2Cum = fsel$AdjR2Cum)
      } else return(NA)
    } else return(NA)
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
  # If we test more than 1 W matrix, 'candidates' only has one class, that is, list. 
  # If we test only one W matrix, 'we have two'candidates' has two classes: listw and nb.
  if (length(class(candidates)) == 1) nbtest <- length(candidates)
  else nbtest <- 1
  
  lenlist <- c()   # Will help with the result output
  
  for (h in 1:k) {
    
    # For model comparison and selection:
    results <- as.data.frame(matrix(nrow = nbtest, ncol = 3))
    colnames(results) <- c("pval_sidak", "R2adj", "NbMEM")     
    # List of the MEM.select matrices (subsets of MEM variables for each signif. W matrix):
    listMEM <- vector("list", nbtest)
    # and corresponding R2:
    listR2 <- vector("list", nbtest)
    # List of 'MEM.test' results:
    listtest <- vector("list", nbtest)
    # List of global W matrices (all MEM variables):
    listW <- vector("list", nbtest)
    
    for (q in 1:nbtest) {
      if (nbtest > 1) W <- scores.listw(candidates[[q]], MEM.autocor = cor[h])
      else W <- scores.listw(candidates, MEM.autocor = cor[h])
      listW[[q]] <- W
      listtest[[q]] <- MEM.test(x, W)
    }
    # Save the results in order to compare them and choose the best model:
    for (i in 1:nbtest) {
      if (is.list(listtest[[i]]) == TRUE) {
        results[i, 1] <- listtest[[i]]$pval
        results[i, 2] <- listtest[[i]]$R2adj
        results[i, 3] <- listtest[[i]]$NbVar
        listMEM[[i]] <- listtest[[i]]$MEM.select
        listR2[[i]] <- listtest[[i]]$AdjR2Cum
      }
    }
    # Selection of the best W matrix (and best eigenvector subset within it):
    if (length(which(results[, 1] <= alpha_thresh)) > 0) {
      best <- which.max(results[, 2])
      if (nbtest == 1) bestlistw = candidates else bestlistw = candidates[[best]]
      if (autocor != "all") {
        lenlist <- c(lenlist, cor[h])
        L <- list(MEM.all = listW[[best]], MEM.select = listMEM[[best]], listw = bestlistw,
                  MEM.AdjR2Cum = listR2[[best]], name = names(candidates)[best], 
                  pval = results[best, 1], R2adj = results[best, 2], 
                  NbVar = results[best, 3], bestw_index = best)
      } else {
        lenlist <- c(lenlist, cor[h])
        if (h == 1) {
          L1 <- list(MEM.all = listW[[best]], MEM.select = listMEM[[best]], 
                     listw = bestlistw, MEM.AdjR2Cum = listR2[[best]], 
                     name = names(candidates)[best], pval = results[best, 1], 
                     R2adj = results[best, 2], NbVar = results[best, 3], bestw_index = best)
        } else {
          L2 <- list(MEM.all = listW[[best]], MEM.select = listMEM[[best]], 
                     listw = bestlistw, MEM.AdjR2Cum = listR2[[best]], 
                     name = names(candidates)[best], pval = results[best, 1], 
                     R2adj = results[best, 2], NbVar = results[best, 3], bestw_index = best)
        }
      }
    }
  }    # End of the for loop
  
  # Output of the MEM.modsel function:
  if (length(lenlist) == 2) {
    cat("\n", "\n", "*****************************************************", "\n",
        "*****************************************************", "\n",
        "A best positive (", sidak, " p-value = ", round(L1$pval, 5), ", R2adj of the", 
        "\n", "selected MEM variables = ", round(L1$R2adj, 3), 
        ") and best negative (", sidak, " p-value = ", round(L2$pval, 5), ",", "\n", 
        "R2adj of the selected MEM variables = ", round(L2$R2adj, 3), ")", 
        "MEM models were selected.", "\n", 
        "*****************************************************", "\n",
        "*****************************************************", "\n",
        "The output of the function is a list of two lists (MEM.pos and MEM.neg).", "\n",
        sep = "")
    list(MEM.pos = L1, MEM.neg = L2)
  } else 
    if (length(lenlist) == 1) { 
      if (autocor == "all") if (lenlist == "positive") L <- L1 else L <- L2
      cat("\n", "\n", "*****************************************************", 
          "\n", "*****************************************************", "\n",
          "A best ", lenlist, " MEM model was selected (", sidak, " p-value = ", 
          round(L$pval, 5), ", R2adj of the selected", "\n",  "MEM variables = ", 
          round(L$R2adj, 3), ").", "\n", 
          "*****************************************************", "\n",
          "*****************************************************", "\n", sep = "")
      list(MEM.all = L$MEM.all, MEM.select = L$MEM.select, listw = L$listw,
           MEM.AdjR2Cum = L$MEM.AdjR2Cum, name = L$name, pval = L$pval, R2adj = L$R2adj, 
           NbVar = L$NbVar, bestw_index = L$bestw_index)
    } else cat("\n", "\n", "*****************************************************", 
               "\n", "*****************************************************", "\n",
               "No significant spatial structure was detected in the data.", "\n",
               "\n", sep = "")   
}
