#' Perform a test of the joint space-environment fraction of a variation partitioning
#' using torus-translation or Moran Spectral Randomisation
#' 
#' This function performs a test of the joint space-environment fraction (JSEF) of a 
#' variation partitioning (VP). The adjusted R-squared (Peres-Neto et al. 2006; R2adj) of 
#' the JSEF is not an actual R2, as it is computed by subtracting the R2adj of other
#' fractions of the VP and therefore has zero degree of freedom. The JSEF
#' can therefore not be computed in the classical way. The function \code{envspace.test}
#' provides two ways of testing this fraction: a torus-translation (or torus randomisation)
#' test (TT), for regular sampling designs, and Moran spectral randomization (MSR). The 
#' function first checks whether the both the environment and the response data display 
#' significant spatial structures, and then proceeds to the significance test if and only 
#' if both datasets display a spatial structure.
#' 
#' @details \code{Ab} can be a vector or a multicolumn matrix or dataframe (multivariate
#' response data). If multivariate, it is greatly advised to transform \code{Ab} prior
#' to performing the VP and testing the JSEF (e.g., Hellinger transformation;
#' Legendre and Gallagher 2001).
#' 
#' \code{MEM_Ab} is a subset of spatial predictors (MEM variables), selected among the
#' complete set of MEM variables of a given spatial weighting matrix (W matrix) chosen
#' by the user. The function \code{\link{listw.candidates}} allows generating a list of 
#' W matrix candidates to be fed to \code{\link{MEM.modsel}}, that optimises the choice 
#' of the W matrix and selects the best subset of spatial predictors within it using the 
#' forward selection with double stopping criterion of Blanchet et al. (2008; see Bauman
#' et al. 2017, and Bauman et al. 2018 in prep. for a review of the methods to select
#' a subset of spatial predictors in spatial eigenvector-based methods, and for a 
#' method of spatial weighting matrix optimisation, respectively).
#' \code{listw_E} corresponds to the spatial weights (W matrix) that will be used to
#' build the MEM variables for the MSR test. The choice of the W matrix for \code{E} can 
#' also be optimised using the \code{\link{MEM.modsel}} function, and attributing the 
#' \code{$listw} output of \code{MEM.modsel} to the argument \code{listw_E} of 
#' \code{envspace.test}.
#' The W matrices selected for \code{Ab} and \code{E) should be chosen separately to 
#' best model the spatial structure of both the response data and the environmental dataset. 
#' If \code{MEM.modsel} is used, as advised, to build \code{MEM_Ab} and \code{listw_E}, 
#' then \code{MEM_Ab} and \code{listw_E} are the \code{$MEM.select} and \code{$listw} parts 
#' of the output of \code{MEM.modsel}, respectively (see example below).
#' 
#' To verify that \code{E} displays a significant spatial pattern, prior to performing the
#' test of the JSEF, a test (see function \code{\link{anova.cca}}) is performed on either 
#' all the MEM variables generated on the basis of \code{listw_E} (\code{autocor_E = "all"), 
#' or on the MEM variables corresponding to a positive or negative autocorrelation 
#' (\code{autocor_E = "positive"} or \code{"negative"}, respectively). If 
#' \code{autocor_E = "all"}, the test is performed separately on the MEM variables 
#' corresponding to positive and negative autocorrelation, and a p-value correction for
#' multiple tests is performed (Sidak correction, see Blanchet et al. 2008). 
#' 
#' \code{E} is a dataset of environmental variables chosen by the user. Although a 
#' significant JSEF may provide an evidence of induced spatial dependence
#' (ISD; Vleminckx et al. 2018), a non significant JSEF only indicate that no
#' induced spatial dependence seems to occur in relation with the chosen environmental
#' variables. This does not exclude the ISD to exist with respect to other unmeasured
#' parameters.
#' 
#' The function needs the environmental variables to be centred and scaled, which is why
#' \code{scale} is set to \code{TRUE} by default. It should only be changed to \code{FALSE}
#' if the user has already scaled \code{E} prior to using \code{envspace.test}.
#' \code{regular} is a logical argument indicating whether a TT test should 
#' be performed additionally to the MSR to test the JSEF. Since, the TT can only
#' be performed on regular sampling designs, \code{regular} should only be set to 
#' \code{TRUE} if the sampling design is a transect or a grid displaying the same number
#' of sites for all lines and columns (although the number of sites can differ between 
#' lines and columns).
#' 
#' @param E Vector, matrix, or dataframe of environmental variables (rows = sites, 
#' columns = variables)
#'            Each variable in E needs to be standardized and normalized prior to the test.
#' @param Ab Vector, matrix, or dataframe of species abundances (rows = sites, 
#' columns = abundances)
#' @param coord Matrix or dataframe of spatial coordinates of the sampled sites
#' @param MEM_Ab Matrix or dataframe of spatial predictors (MEM variables) for the response
#' data (\code{Ab})
#' @param listw_E An object of class \code{listw} (spatial weights) created by the functions
#' of the \bold{\code{spdep}} package or returned by the function \code{\link{MEM.modsel}} in
#' its \code{$listw} output
#' @param autocor_E A string indicating whether all the MEM variables built on the basis of
#' \code{listW_E} should be used to test the significance of a spatial structure in \code{E},
#' (\code{"all"}), or only those corresponding to positive (\code{"positive"}) or 
#' negative (\code{"negative"}) autocorrelation; Default is "positive" 
#' @param scale Logical value indicating whether the environmental variables should be
#' centered and scaled (standardised data are necessary); Default is \code{TRUE}
#' @param regular Logical argument indicating whether a torus-translation test is
#' performed, additionnaly to the MSR. Set to \code{TRUE} if the sampling design is regular 
#' (same number of sites on each line, same number of sites on each column), set to
#' \code{FALSE} otherwise; Default is \code{FALSE}
#' @param ntoro Number of permutations performed using torus-translations; Default is 999
#' @param nMSR  Number of permutations performed using MSR; Default is 999
#' @param MSRmethod Algorithm of \code{\link{msr}} to be used to perform the MSR. The
#' three available procedures are "singleton" (default), "pair", and "triplet" (see
#' \code{\link{msr}} for details)
#' @param alpha Threshold value of null hypothesis rejection for the test of a 
#' spatial structure in the environment, and for the joint environment-space fraction of 
#' the VP; Default is 0.05
#' @param summary Logical; Whether a message summarising the results should be returned
#' or not; Default is \code{TRUE}
#' 
#' @return The function returns a list containing the following elements: \describe{
#' \item{R2adj}{The adjusted R-squared value of the JSEF.}
#' \item{pval_TT}{The significance value of the JSEF obtained by TT.}
#' \item{pval_MSR}{The significance value of the JSEF obtained by MSR.}
#' }
#' 
#' @author Jason Vleminckx and David Bauman, \email{jasvlx86@@gmail.com}, 
#' \email{davbauman@@gmail.com}
#' 
#' @seealso \code{\link{varpart}}, \code{\link{MEM.modsel}}, \code{\link{listw.candidates}}
#' 
#' @references 
#' Bauman D., Drouet T., Dray S. and Vleminckx J. (2017) Disentangling good from bad 
#' practices in the selection of spatial or phylogenetic eigenvectors. Ecography 
#' 
#' Blanchet G., Legendre P. and Borcard D. (2008) Forward selection of explanatory variables.
#' Ecology, 89(9), 2623--2632
#' 
#' Legendre P., Gallagher E.D. (2001) Ecologically meaningful transformations for 
#' ordination of species data. Oecologia, 129(2), 271--280
#' 
#' Peres-Neto P., Legendre P., Dray S., Borcard D. (2006) Variation partitioning of 
#' species data matrices: estimation and comparison of fractions. Ecology, 87(10), 
#' 2614--2625
#' 
#' Peres-Neto P. and Legendre P. (2010) Estimating and controlling for spatial structure 
#' in the study of ecological communities. Global Ecology and Biogeography, 19, 174--184
#' 
#' @keywords spatial
#' 
#' @examples
#' if(require(vegan)) { 
#' # Illustration of the test of the JSEF on the oribatid mite data
#' # (Borcard et al. 1992, 1994 for details on the dataset):
#' # Community data (response matrix):
#' data(mite)
#' # Hellinger-transformation of the community data (Legendre and Gallagher 2001):
#' Y <- decostand(mite, method = "hellinger")
#' # Environmental explanatory dataset:
#' data(mite.env)
#' # For the example, we only use two numerical explanatory variables:
#' env <- mite.env[, 1:2]
#' dim(Y)
#' dim(env)
#' # Coordinates of the 70 sites:
#' data(mite.xy)
#' coord <- mite.xy
#' 
#' ### Building a list of spatial weighting matrix (W matrix) candidates for the 
#' ### optimisation of the W matrix selection for 'Y' and 'env':
#' # We create five candidates: a connectivity matrix based on a Gabriel graphs, on
#' # a minimum spanning tree (i.e., two contrasted graph-based W matrices), either
#' # not weighted, or weighted by a linear function decreasing with the distance),
#' # and a distance-based W matrix corresponding to the connectivity and weighting
#' # criteria of the original PCNM method:
#' candidates <- listw.candidates(coord, del = FALSE, gab = TRUE, rel = FALSE,
#'                                mst = TRUE, DB = FALSE, PCNM = TRUE, 
#'                                binary = TRUE, flin = TRUE, fconcdown = FALSE,
#'                                fconcup = FALSE)
#' ### Optimisation of the choice of a W matrix:
#' # W matrix for 'Y':
#' modsel_Y <- MEM.modsel(Y, candidates, correction = TRUE, autocor = "positive")
#' 
#' paste("The best W matrix for 'Y' was ", modsel_Y$name, ". The forward selection ",
#'       "with double stopping criterion (Blanchet et al. 2008) selected a best subset",
#'       " of spatial predictors within this W matrix. This subset contains ", 
#'       modsel_Y$NbVar, " MEM variables and has an adjusted R2 of ", 
#'       round(modsel_Y$R2adj, 3), ".", sep = "")
#' 
#' modsel_env <- MEM.modsel(env, candidates, correction = TRUE, autocor = "all")
#' 
#' paste("The best W matrix for 'env' was ", modsel_env$name, ". The forward selection ",
#'       "with double stopping criterion (Blanchet et al. 2008) selected a best subset",
#'       " of spatial predictors within this W matrix. This subset contains ", 
#'       modsel_env$NbVar, " MEM variables and has an adjusted R2 of ", 
#'       round(modsel_env$R2adj, 3), ".", sep = "")
#' 
#' ### We perform the variation partitioning:
#' # Subset of selected MEM variables within the best W matrix:
#' MEM_Ab <- modsel_Y$MEM.select
#' 
#' VP <- varpart(Y, env, MEM_Ab)
#' plot(VP)
#' 
#' # Test of the joint space-environment fraction (fraction [b]):
#' JSEF.test <- envspace.test(E = env, Ab = Y, coord = coord, MEM_Ab = MEM_Ab,
#'                            listw_E = modsel_env$listw, scale = TRUE, regular = FALSE,
#'                            nMSR = 999)
#' JSEF.test$R2adj
#' JSEF.test$pval_MSR
#' 
#' # The JSEF is very highly significant (p-value = 0).
#' }
#' 
#' @importFrom vegan rda anova.cca RsquareAdj varpart
#' @export

"envspace.test" <- function(E, 
                            Ab, 
                            coord, 
                            MEM_Ab,
                            listw_E,
                            autocor_E = c("positive", "negative", "all"),
                            scale = TRUE, 
                            regular = FALSE, 
                            ntoro = 999, 
                            nMSR = 999,
                            MSRmethod = "singleton",
                            alpha = 0.05,
                            summary = TRUE) {
  
  Ab <- as.data.frame(Ab)
  coord <- as.matrix(coord)
  
  if (any(is.na(Ab)) | any(is.na(coord)))
    stop("NA entries in x or coord")
  if (nrow(Ab) != nrow(coord))
    stop("different number of rows")
  if (inherits(listw_E, what = "listw") == FALSE)
    stop("listw_E is not an object of class 'listw'")
  
  if(is.vector(E) == "TRUE") 
    nE <- 1
    else nE <- ncol(E)
  if(is.vector(Ab) == "TRUE") 
    nAb <- 1
    else nAb <- ncol(Ab)
    
  M <- as.data.frame(matrix(0, ncol = (ncol(coord) + nAb + nE), nrow = nrow(cbind(Ab, E))))
  
  M[, c(1: 2)]               <- coord
  M[, c(3:(2 + nAb))]        <- Ab
  M[, c((3 + nAb): ncol(M))] <- E

  # Testing whether the environment displays a significant spatial structure:
  # *************************************************************************
  autocor_E <- match.arg(autocor_E)
  
  if (autocor_E != "all") {
    MEM_E <- scores.listw(listw_E, MEM.autocor = autocor_E)
    test.envspa <- anova.cca(rda(E, MEM_E))
    pval <- test.envspa$Pr[1]
  } else {
    pval <- c()
    for (i in 1:2) {
      if (i == 1) aut <- "positive" else aut <- "negative"
      MEM_E <- scores.listw(listw_E, MEM.autocor = aut)
      test.envspa <- anova.cca(rda(E, MEM_E))
      p <- 1-(1-test.envspa$Pr[1])^2   # Sidak correction
      pval <- c(pval, p)
    }
  }
  
  if (length(which(pval <= alpha)) == 0) {
    return("No spatial structure detected in the environment;",
           "The joint environment-space fraction should not be considered.")
  }
  
  # Testing whether the response data is significantly related to the environment:
  # ******************************************************************************

  test.env <- anova.cca(rda(Ab, E))
  if(test.env$Pr[1] > alpha) {
    return("No significant relation between the response data and the environment;", 
           "The joint environment-space fraction should not be considered.")
  }
  
  #*********************************************************************************
  ###   TESTING THE R2 of the env-space FRACTION OF THE VARIATION PARTITIONING   ###
  #*********************************************************************************
  
  ## R2env-space value:
  #********************
  R2.b <- varpart(M[, c(3:(2+nAb))], E, MEM_Ab)$part$indfract$Adj.R.square[2]

     #***********************************
     ###   USING Torus-translations   ###
     #***********************************
  
  if (regular == "TRUE") {
    
    ## vector to be filled with null values of R2env-space:
    #******************************************************
    E.b.T <- c()
    ## Calculate p-value:
    #********************
    for (k in 1:ntoro) {
      M2 <- M
      sens <- runif(1, min = 0, max = 1)
      rx <- ceiling(runif(1, 0.01, 0.99) * 50)
      ry <- ceiling(runif(1, 0.01, 0.99) * 25)
      if (sens <= 0.5) { 
        M2[, 1] <- (M2 [, 1] + rx) %% max(M2 [, 1])
        M2[, 2] <- (M2 [, 2] + ry) %% max(M2 [, 2])
      } else { 
        M2[, 1] <- max(M2[, 1]) - (M2[, 1] + rx) %% max(M2[, 1])
        M2[, 2] <- max(M2[, 2]) - (M2[, 2] + ry) %% max(M2[, 2]) 
      }
      M3 <- M
      M2 <- M2[order(M2[, 2]), ]
      M2 <- M2[order(M2[, 1]), ]
      M3[, c((3+nAb): ncol(M))] <- M2[, c((3+nAb): ncol(M))]
      E.toro <- M3[, c((3+nAb):(3+nE))]
      E.b.T <- c(E.b.T, varpart(M[, c(3:(2+nAb))], E.toro, 
                                MEM_Ab)$part$indfract$Adj.R.square[2])
    }  # end of the "for (k in 1:ntoro) {"
    
    p.value.toro <- length(E.b.T[E.b.T > R2.b]) / ntoro
    
  } else p.value.toro <- "Not computed due to inappropriate sampling design"
  
  #******************************************************************
  ###   TESTING R2env-space USING MORAN SPECTRAL RANDOMISATIONS   ###
  #******************************************************************

  ## Moran Spectral Randomisation of environmental structures (N = nMSR):
  MSR.ENV <- msr(E, listw_E, nrepet = nMSR, method = MSRmethod) 
  
  ## vector to be filled with null values of R?env-space:
  #******************************************************
  E.b.MSR <- c()
  
  ## calculate p-value:
  #********************
  for(k in 1:nMSR){
      E.b.MSR <- c(E.b.MSR, varpart(M[, c(3:(2+nAb))], 
                                    MSR.ENV[[k]], MEM_Ab)$part$indfract$Adj.R.square[2])
  }
  
  p.value.MSR <- length(E.b.MSR[E.b.MSR > R2.b]) / nMSR

  if (summary == "TRUE") {
    if (regular == TRUE)
      cat("\n", "************************************************************************", 
          "\n", "The R2 value of the joint space-environment fraction (R2adj) ", 
          "and the corresponding", "\n", "p-values computed with the TT test and the",
          " MSR are:", "\n", "R2adj = ", round(R2.b, 3), ", p-value TT = ", p.value.toro, 
          ", and p-value MSR = ", p.value.MSR, ".",
          "\n", "************************************************************************", 
          "\n", "\n", sep = "")
    else
      cat("\n", "************************************************************************", 
          "\n", "The R2 value of the joint space-environment fraction (R2adj) ", 
          "and the corresponding", "\n", "p-value computed with the MSR test are:",
          "\n", "R2adj = ", round(R2.b, 3), ", p-value MSR = ", p.value.MSR, ".", 
          "\n", "************************************************************************", 
          "\n", "\n", sep = "")
  }
  
  list(R2adj = R2.b, pval_TT = p.value.toro, pval_MSR = p.value.MSR)
}
