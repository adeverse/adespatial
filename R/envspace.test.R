#' Perform a test of the shared space-environment fraction of a variation partitioning
#' using torus-translation (TT) or Moran Spectral Randomisation (MSR)
#' 
#' The function uses two different spatially-constrained null models to test the shared 
#' space-environment fraction (SSEF, or fraction [b]) of a variation partitioning of two 
#' explanatory components.
#' 
#' @details The function tests the SSEF (also known as fraction [b]) of a variation 
#' partitioning of a response variable or matrix (\code{y}) between an environmental and a
#' spatial component (\code{env}, and \code{MEM.spe}, respectively). The SSEF is the 
#' explained variation of \code{y} shared by \code{env} and \code{MEM.spe}. 
#' The adjusted R-squared (Peres-Neto et al. 2006; R2adj) of the SSEF is not an 
#' actual R2, as it is computed by subtracting the adjusted R2adj of other fractions and 
#' therefore has zero degree of freedom (Legendre and Legendre 2012). 
#' The SSEF can therefore not be computed in the classical way (residuals permutation; 
#' Anderson and Legendre 1999, Legendre and Legendre 2012). 
#' 
#' The function \code{envspace.test} provides two ways of testing this fraction, that is,
#' spatially-constrained null models based either on a torus-translation test (TT) (for 
#' regular sampling designs only), or on Moran spectral randomizations (MSR) (for any type
#' of sampling design). The test of the SSEF should only be performed if both the global 
#' models of \code{y} against all the environmental variables and against all spatial variables
#' are significant (see Bauman et al. 2018c).
#' The function first checks whether the environment displays significant spatial structures, 
#' and then proceeds to the test of the SSEF if this condition is fulfilled (details in
#' Bauman et al. 2018c).
#' 
#' \code{spe} can be a vector or a multicolumn matrix or dataframe (multivariate
#' response data). If multivariate, it is greatly advised to transform \code{spe} prior
#' to performing the variation partitioning and testing the SSEF (e.g., Hellinger 
#' transformation; see Legendre and Gallagher 2001).
#' 
#' \code{MEM.spe} is a set of spatial predictors (MEM variables). It is recommended to be
#' a well-defined subset of MEM variables selected among the complete set generated from 
#' the spatial weighting matrix (SWM) (see review about spatial eigenvector selection in
#' Bauman et al. 2018a). 
#' Optimising the selection of a subset of forward-selected MEM variables
#' among a set of candidate SWMs has been shown to increase statistical power as well as
#' R2-estimation accuracy (Bauman et al. 2018b). To do so, \code{MEM.spe} can be generated
#' using \code{\link{listw.candidates}} followed by \code{\link{listw.select}}. If a SWM has
#' already been selected in another way, then \code{\link{mem.select}} can be used to
#' generate the MEM variables and to select an optimal subset among them, which can then 
#' be used as \code{MEM.spe} in \code{envspace.test} (see \code{Details} of function 
#' \code{mem.select}).
#' \code{listw.env} corresponds to the SWM that will be used to test for a spatial structure
#' in \code{env}, and to build the MEM variables for the MSR test. 
#' The choice of the SWM for \code{env} can also be optimised with \code{\link{listw.select}}.
#' The SWMs selected for \code{spe} and \code{env} should be optimised separately to 
#' best model the spatial structure of both \code{spe} and \code{env} (see example).
#' 
#' To verify that \code{env} displays a significant spatial pattern, prior to performing the
#' test of the SSEF, a residuals permutation test is performed on the global set of MEM 
#' variables (generated internally from \code{listw.env}) associated to the type of 
#' spatial structure of interest (see argument \code{MEM.autocor}). This test is performed 
#' with \code{mem.select}. The choice of \code{MEM.autocor} should be made according to 
#' the \code{MEM.autocor} argument used to build \code{MEM.spe}. 
#' 
#' \code{env} is a dataset of environmental variables chosen by the user. We recommend dealing
#' with collinearity issues prior to performing the variation partitioning and the test of
#' the SSEF (see Dormann et al. 2013 for a review of methods to cope with collinearity).
#' 
#' The function needs the environmental variables to be centred and scaled, which is why
#' \code{scale} is set to \code{TRUE} by default. It should only be changed to \code{FALSE}
#' if the user has already scaled \code{env} prior to using \code{envspace.test}.
#' \code{regular} is a logical argument indicating whether a TT test should 
#' be performed additionally to the MSR to test the SSEF. Since the TT can only
#' be performed on regular sampling designs, \code{regular} should only be set to 
#' \code{TRUE} if the sampling design is either a transect, or a grid displaying the 
#' same number of sites for all lines and columns (although the number of sites per column 
#' can differ from the number of sites per line).
#' 
#' \code{listw.env} is the SWM used by the MSR to generate spatially-constrained null
#' environmental variables. It should ideally be a SWM optimised on the basis of \code{env}  
#' using the function \code{listw.select}, with the argument \code{method = "global"} (see
#' \code{Details} of function \code{mem.select} for an explanation). 
#' This will allow detecting the spatial structures of \code{env} as accurately as possible, 
#' hence allowing MSR to generate null environmental variables as spatially faithful to the 
#' original ones. 
#' It is also on the basis of \code{listw.env} that MEM variables will be generated to test
#' whether \code{env} is spatially structured (i.e. global test) prior to perform the test of
#' the SSEF.
#' 
#' It is worth mentioning that, although a significant SSEF may provide evidence of an 
#' induced spatial dependence (Bauman et al. 2018c), a non-significant SSEF only indicates 
#' that no induced spatial dependence could be detected in relation with the chosen 
#' environmental variables. This does not exclude that this effect may exist with respect 
#' to some unmeasured variables.
#' 
#' @param env Vector, matrix, or dataframe of environmental variables (rows = sites, 
#' columns = variables)
#' @param spe Vector, matrix, or dataframe of response variable(s) (e.g. species abundances)
#' @param coord Matrix or dataframe of spatial coordinates of the sampled sites
#' @param MEM.spe Matrix or dataframe of spatial predictors (MEM variables) selected for 
#' \code{spe}
#' @param listw.env An object of class \code{listw} (spatial weights) created by the functions
#' of the \code{spdep} package or returned by \code{\link{listw.candidates}}
#' @param MEM.autocor A string indicating the type of spatial structure of interest for 
#' \code{env} (\code{"positive"}, \code{"negative"}, or \code{"all"}, for positive, negative, 
#' or both types of spatial autocorrelations, respectively); Default is \code{"positive"}
#' @param scale Logical value indicating whether the environmental variables should be
#' centered and scaled (standardised data are necessary); Default is \code{TRUE}
#' @param regular Logical argument indicating whether a torus-translation test will be
#' performed, in addition to the MSR. Set to \code{TRUE} only if the sampling design is regular 
#' (same number of sites on each line, same number of sites on each column). Set to
#' \code{FALSE} otherwise; Default is \code{FALSE}
#' @param nperm  Number of permutations performed; Default is 999
#' @param MSR.method Algorithm of \code{\link{msr}} to be used to perform the MSR. The
#' three available procedures are \code{"singleton"} (default), \code{"pair"}, and 
#' \code{"triplet"} (see \code{\link{msr}} for details)
#' @param alpha Threshold value of null hypothesis rejection for the test of a 
#' spatial structure in the environment, and for the shared environment-space fraction of 
#' the variation partitioning; Default is 0.05
#' @param summary Logical; Whether a message summarising the results should be returned
#' or not; Default is \code{FALSE}
#' 
#' @return If the condition of \code{env} being spatially structured is fulfilled, the test 
#' is performed and the function returns a list containing the following elements: 
#' \describe{
#' \item{R2adj}{The adjusted R-squared value of the SSEF.}
#' \item{pval_TT}{The significance value of the SSEF obtained by TT.}
#' \item{pval_MSR}{The significance value of the SSEF obtained by MSR.}
#' }
#' Otherwise, the function returns a message informing why the test was not performed.
#' 
#' @author David Bauman and Jason Vleminckx, \email{davbauman@@gmail.com}, 
#' \email{jasvlx86@@gmail.com}
#' 
#' @seealso \code{\link{varpart}}, \code{\link{listw.select}}, \code{\link{listw.candidates}}, \code{\link{mem.select}}
#' 
#' @references Anderson M. and Legendre P. (1999) An empirical comparison of permutation 
#' methods for tests of partial regression coefficients in a linear model. Journal of
#' Statistical Computation and Simulation, 62(3), 271--303
#' 
#' Bauman D., Drouet T., Dray S. and Vleminckx J. (2018a) Disentangling good 
#' from bad practices in the selection of spatial or phylogenetic eigenvectors. Ecography, 
#' 41, 1--12
#' 
#' Bauman D., Fortin M-J, Drouet T. and Dray S. (2018b) Optimizing the choice 
#' of a spatial weighting matrix in eigenvector-based methods. Ecology
#' 
#' Bauman D., Vleminckx J., Hardy O., Drouet T. (2018c) Testing and interpreting the 
#' shared space-environment fraction in variation partitioning analyses of ecological data.
#' Oikos
#' 
#' Blanchet G., Legendre P. and Borcard D. (2008) Forward selection of explanatory variables.
#' Ecology, 89(9), 2623--2632
#' 
#' Legendre P., Gallagher E.D. (2001) Ecologically meaningful transformations for 
#' ordination of species data. Oecologia, 129(2), 271--280
#' 
#' Legendre P. and Legendre L. (2012) Numerical Ecology, Elsevier, Amsterdam
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
#' # Illustration of the test of the SSEF on the oribatid mite data
#' # (Borcard et al. 1992, 1994 for details on the dataset):
#' # Community data (response matrix):
#' data(mite)
#' # Hellinger-transformation of the community data (Legendre and Gallagher 2001):
#' Y <- decostand(mite, method = "hellinger")
#' # Environmental explanatory dataset:
#' data(mite.env)
#' # We only use two numerical explanatory variables:
#' env <- mite.env[, 1:2]
#' dim(Y)
#' dim(env)
#' # Coordinates of the 70 sites:
#' data(mite.xy)
#' coord <- mite.xy
#' 
#' ### Building a list of candidate spatial weighting matrices (SWMs) for the 
#' ### optimisation of the SWM selection, separately for 'Y' and 'env':
#' # We create five candidate SWMs: a connectivity matrix based on a Gabriel graphs, on
#' # a minimum spanning tree (i.e., two contrasted graph-based SWMs), either
#' # not weighted, or weighted by a linear function decreasing with the distance),
#' # and a distance-based SWM corresponding to the connectivity and weighting
#' # criteria of the original PCNM method:
#' candidates <- listw.candidates(coord, nb = c("gab", "mst", "pcnm"), weights = c("binary",
#'                                                                                 "flin"))
#' ### Optimisation of the selection of a SWM:
#' # SWM for 'Y' (based on the best forward-selected subset of MEM variables):
#' modsel.Y <- listw.select(Y, candidates, method = "forward", MEM.autocor = "positive",
#'                          p.adjust = TRUE)
#'                          
#' names(candidates)[modsel.Y$best.id]                 # Best SWM selected
#' modsel.Y$candidates$Pvalue[modsel.Y$best.id]        # Adjusted p-value of the global model
#' modsel.Y$candidates$N.var[modsel.Y$best.id]         # Nb of forward-selected MEM variables
#' modsel.Y$candidates$R2Adj.select[modsel.Y$best.id]  # Adjusted R2 of the selected MEM var.
#' 
#' # SWM for 'env' (method = "global" for the optimisation, as all MEM variables are required
#' # to use MSR):
#' modsel.env <- listw.select(env, candidates, method = "global", MEM.autocor = "positive",
#'                            p.adjust = TRUE)
#'
#' names(candidates)[modsel.env$best.id]                  # Best SWM selected
#' modsel.env$candidates$Pvalue[modsel.env$best.id]       # Adjusted p-value of the global model
#' modsel.env$candidates$N.var[modsel.env$best.id]        # Nb of forward-selected MEM variables
#' modsel.env$candidates$R2Adj.select[modsel.env$best.id] # Adjusted R2 of the selected MEM var.
#' 
#' ### We perform the variation partitioning:
#' # Subset of selected MEM variables within the best SWM:
#' MEM.spe <- modsel.Y$best$MEM.select
#' 
#' VP <- varpart(Y, env, MEM.spe)
#' plot(VP)
#' 
#' # Test of the shared space-environment fraction (fraction [b]):
#' SSEF.test <- envspace.test(spe, env, coord, MEM.spe, 
#'                            listw.env = candidates[[modsel.env$best.id]], scale = TRUE, 
#'                            regular = FALSE, nperm = 999)
#' SSEF.test$R2adj
#' SSEF.test$pval_MSR
#' 
#' # The SSEF is highly significant, indicating a potential induced spatial dependence.
#' }
#' 
#' @importFrom vegan rda anova.cca RsquareAdj varpart
#' @export

"envspace.test" <- function(spe, 
                            env, 
                            coord, 
                            MEM.spe,
                            listw.env,
                            MEM.autocor = c("positive", "negative", "all"),
                            scale = TRUE, 
                            regular = FALSE, 
                            nperm = 999,
                            MSR.method = "singleton",
                            alpha = 0.05,
                            summary = FALSE) {
  
  spe <- as.data.frame(spe)
  coord <- as.matrix(coord)
  
  if (any(is.na(spe)) | any(is.na(coord)))
    stop("NA entries in spe or coord")
  if (nrow(spe) != nrow(coord))
    stop("different number of rows")
  if (inherits(listw.env, what = "listw") == FALSE)
    stop("listw.env is not an object of class 'listw'")
  if (is.vector(MEM.spe) == FALSE) {
    if (ncol(MEM.spe) == nrow(coord)-1)
      stop("n-1 MEM variables in 'MEM.spe'. Select a subset to avoid saturated model issues.")
  }

  nspe <- ncol(spe)  
  if(is.vector(env) == "TRUE") 
    nenv <- 1
  else nenv <- ncol(env)
    
  M <- as.data.frame(matrix(0, ncol = (ncol(coord) + nspe + nenv), 
                            nrow = nrow(cbind(spe, env))))
  
  M[, c(1: 2)] <- coord
  M[, c(3:(2 + nspe))] <- spe
  M[, c((3 + nspe): ncol(M))] <- env

  # Testing whether the environment displays a significant spatial structure:
  # *************************************************************************
  MEM.autocor <- match.arg(MEM.autocor)
  global <- mem.select(env, listw.env, method = "global", 
                       MEM.autocor = MEM.autocor, alpha = alpha)$global.test
  if(MEM.autocor != "all") 
    pval <- global$pvalue 
  else 
    pval <- c(global$positive$pvalue, global$negative$pvalue)
  
  if (length(which(pval <= alpha)) == 0) {
    message("No significant spatial structure detected in 'env'; the test was not performed")
    return()
  }
  
  # Testing whether the response data displays a significant spatial structure:
  # ***************************************************************************
  test.spe <- anova.cca(rda(spe, MEM.spe))
  if(test.spe$Pr[1] > alpha) {
    message("No significant spatial structure detected in 'spe'; the test was not performed")
    return()
  }
  
  # Testing whether the response data is significantly related to the environment:
  # ******************************************************************************
  test.env <- anova.cca(rda(spe, env))
  if(test.env$Pr[1] > alpha) {
    message("No significant relation between 'spe' and 'env'; the test was not performed")
    return()
  }
  
  #*********************************************************************************
  ###   TESTING THE R2 of the env-space FRACTION OF THE VARIATION PARTITIONING   ###
  #*********************************************************************************
  
  ## Define whether the SSEF computed from the unadjusted fractions is negative:
  SSEF.unadj <- RsquareAdj(rda(spe, env))$r.squared - RsquareAdj(rda(spe, env, MEM.spe))$r.squared
  alternative <- ifelse(SSEF.unadj < 0, "smaller", "greater")
  
  ## Observed SSEF value:
  #**********************
  R2.b <- varpart(M[, c(3:(2+nspe))], env, MEM.spe)$part$indfract$Adj.R.square[2]

     #***********************************
     ###   USING Torus-translations   ###
     #***********************************
  
  if (regular == "TRUE") {
    
    ## vector to be filled with null values of R2env-space:
    #******************************************************
    E.b.T <- c()
    ## Calculate p-value:
    #********************
    for (k in 1:nperm) {
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
      M3[, c((3 + nspe): ncol(M))] <- M2[, c((3 + nspe): ncol(M))]
      E.toro <- M3[, c((3 + nspe):(3 + nenv))]
      E.b.T <- c(E.b.T, varpart(M[, c(3:(2+nspe))], E.toro, 
                                MEM.spe)$part$indfract$Adj.R.square[2])
    }  # end of the "for (k in 1:nperm) {"
    
    if (alternative == "greater")
      p.value.toro <- length(which(c(R2.b, E.b.T) >= R2.b)) / (nperm+1)
    else p.value.toro <- length(which(c(R2.b, E.b.T) <= R2.b)) / (nperm+1)
    
  } else p.value.toro <- "Not computed"
  
  #******************************************************************
  ###   TESTING R2env-space USING MORAN SPECTRAL RANDOMISATIONS   ###
  #******************************************************************

  ## Moran Spectral Randomisation of environmental structures:
  MSR.ENV <- msr(env, listw.env, nrepet = nperm, method = MSR.method, simplify = FALSE) 
  
  ## vector to be filled with null values of R2env-space:
  #******************************************************
  E.b.MSR <- c()
  
  ## calculate p-value:
  #********************
  for(k in 1:nperm){
      E.b.MSR <- c(E.b.MSR, varpart(M[, c(3:(2+nspe))], 
                                    MSR.ENV[[k]], MEM.spe)$part$indfract$Adj.R.square[2])
  }
  
  if (alternative == "greater")
    p.value.MSR <- length(which(c(R2.b, E.b.MSR) >= R2.b)) / (nperm+1)
  else p.value.MSR <- length(which(c(R2.b, E.b.MSR) <= R2.b)) / (nperm+1)

  if (summary == "TRUE") {
    if (regular == TRUE)
      cat("\n", "************************************************************************", 
          "\n", "The R2 value of the shared space-environment fraction (R2adj) ", 
          "and the corresponding", "\n", "p-values computed with the TT test and the",
          " MSR are:", "\n", "R2adj = ", round(R2.b, 3), ", p-value TT = ", p.value.toro, 
          ", and p-value MSR = ", p.value.MSR, ".",
          "\n", "************************************************************************", 
          "\n", "\n", sep = "")
    else
      cat("\n", "************************************************************************", 
          "\n", "The R2 value of the shared space-environment fraction (R2adj) ", 
          "and the corresponding", "\n", "p-value computed with the MSR test are:",
          "\n", "R2adj = ", round(R2.b, 3), ", p-value MSR = ", p.value.MSR, ".", 
          "\n", "************************************************************************", 
          "\n", "\n", sep = "")
  }
  
  list(R2adj = R2.b, pval_TT = p.value.toro, pval_MSR = p.value.MSR)
}
