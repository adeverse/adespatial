#' Perform a test of the shared space-environment fraction of a variation partitioning
#' using torus-translation (TT) or Moran Spectral Randomisation (MSR)
#'
#' The function uses two different spatially-constrained null models to test the shared
#' space-environment fraction (SSEF) of a variation partitioning of two explanatory
#' components. In the notation of Borcard et al. (1992) and Legendre and Legendre (2012),
#' the SSEF is fraction [b]; note that \code{\link[vegan]{varpart}} labels this shared
#' fraction [c] since vegan 2.6-4 (see \code{Details}).
#'
#' @details The function tests the SSEF of a variation partitioning of a response variable
#' or matrix (\code{spe}) between an environmental and a spatial component (\code{env}, and
#' \code{MEM.spe}, respectively). The SSEF is the explained variation of \code{spe} shared
#' by \code{env} and \code{MEM.spe}.
#' The adjusted R-squared (Peres-Neto et al. 2006; R2adj) of the SSEF is not an
#' actual R2, as it is computed by subtracting the adjusted R2adj of other fractions and
#' therefore has zero degree of freedom (Legendre and Legendre 2012).
#' The SSEF can therefore not be computed in the classical way (residuals permutation;
#' Anderson and Legendre 1999, Legendre and Legendre 2012).
#'
#' The SSEF is computed directly from its definition, as
#' R2adj(\code{env}) + R2adj(\code{MEM.spe}) - R2adj(\code{env} + \code{MEM.spe}),
#' where each R2adj is the adjusted R-squared of a redundancy analysis (RDA) of \code{spe}
#' (\code{\link[vegan]{rda}} and \code{\link[vegan]{RsquareAdj}}). This value is identical
#' to the shared fraction returned by \code{\link[vegan]{varpart}}, but it does not depend
#' on the order or labelling of the fractions in the output of \code{varpart}.
#' These changed in vegan 2.6-4 (October 2022): with two explanatory tables, the
#' shared fraction, historically labelled [b] and reported in second position (Borcard et
#' al. 1992), is now labelled [c] and reported in third position, after the two unique
#' fractions (now [a] and [b]).
#'
#' The function \code{envspace.test} provides two ways of testing this fraction, that is,
#' spatially-constrained null models based either on a torus-translation test (TT) (for
#' regular sampling designs only), or on Moran spectral randomizations (MSR) (for any type
#' of sampling design). In both cases, the response (\code{spe}) and the spatial predictors
#' (\code{MEM.spe}) are kept fixed, and null environmental datasets that preserve the
#' spatial structure of \code{env} are generated to build the null distribution of the
#' SSEF. The test of the SSEF should only be performed if both the global
#' models of \code{spe} against all the environmental variables and against all spatial
#' variables are significant (see Bauman et al. 2018c).
#' The function first checks whether the environment displays significant spatial
#' structures, and then proceeds to the test of the SSEF if this condition is fulfilled
#' (details in Bauman et al. 2018c).
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
#' \code{env} is a dataset of numeric environmental variables chosen by the user. We
#' recommend dealing with collinearity issues prior to performing the variation
#' partitioning and the test of the SSEF (see Dormann et al. 2013 for a review of methods
#' to cope with collinearity).
#'
#' \code{regular} is a logical argument indicating whether a TT test should
#' be performed instead of the MSR to test the SSEF. Since the TT can only
#' be performed on regular sampling designs, \code{regular} should only be set to
#' \code{TRUE} if the sampling design is either a transect, or a complete grid of evenly
#' spaced sites displaying the same number of sites for all lines and columns (although the
#' number of sites per column can differ from the number of sites per line). The grid
#' positions of the sites are derived from \code{coord}, so that the rows of the data can be
#' in any order and the coordinates can be expressed in any unit. Each torus translation
#' shifts the environmental data by a random number of grid cells drawn uniformly over the
#' whole extent of each axis (with wrap-around), and is combined with a point reflection
#' (180-degree rotation) of the grid with probability 0.5. The function stops with an error
#' if \code{coord} does not describe such a regular design.
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
#' @note Up to adespatial 0.3-29, \code{envspace.test} extracted the SSEF from the output
#' of \code{\link[vegan]{varpart}} by position (second row of \code{$part$indfract}).
#' With vegan >= 2.6-4, this row is the pure spatial fraction (\code{MEM.spe} | \code{env}),
#' so that these versions tested the pure spatial fraction instead of the SSEF. Because
#' R2adj(\code{MEM.spe}) is constant across permutations, the pure spatial fraction is a
#' decreasing linear function of the SSEF, and the erroneous test was the exact mirror of
#' the correct one: the reported observed value was the pure spatial fraction, the sign of
#' the standardized observation was inverted, and, in the absence of ties, the correct
#' p-value is p = (nperm + 2) / (nperm + 1) - p.reported. This applies to the MSR test
#' and to the TT test with a univariate response. The TT test of these versions was also
#' incorrect for multivariate responses, regardless of the vegan version (the null
#' environmental datasets included columns of the response data); in addition, its torus
#' translations were only valid when the rows were ordered by x then y coordinates and the
#' coordinates did not start at 0, and the shifts were limited to hard-coded maxima of 50
#' (x) and 25 (y) coordinate units.
#'
#' @param spe Vector, matrix, or dataframe of response variable(s) (e.g. species abundances)
#' @param env Vector, matrix, or dataframe of numeric environmental variables (rows = sites,
#' columns = variables)
#' @param coord Matrix or dataframe of spatial coordinates of the sampled sites
#' @param MEM.spe Matrix or dataframe of spatial predictors (MEM variables) selected for
#' \code{spe}
#' @param listw.env An object of class \code{listw} (spatial weights) created by the functions
#' of the \code{spdep} package or returned by \code{\link{listw.candidates}}
#' @param MEM.autocor A string indicating the type of spatial structure of interest for
#' \code{env} (\code{"positive"}, \code{"negative"}, or \code{"all"}, for positive, negative,
#' or both types of spatial autocorrelations, respectively); Default is \code{"positive"}
#' @param regular Logical argument indicating whether a torus-translation test will be
#' performed instead of the MSR. Set to \code{TRUE} only if the sampling design is regular
#' (complete grid of evenly spaced sites, or transect). Set to \code{FALSE} otherwise;
#' Default is \code{FALSE}
#' @param nperm  Number of permutations performed; Default is 999
#' @param MSR.method Algorithm of \code{\link{msr}} to be used to perform the MSR. The
#' three available procedures are \code{"singleton"} (default), \code{"pair"}, and
#' \code{"triplet"} (see \code{\link{msr}} for details)
#' @param alpha Threshold value of null hypothesis rejection for the test of a
#' spatial structure in the environment, and for the shared environment-space fraction of
#' the variation partitioning; Default is 0.05
#'
#' @return If the condition of \code{env} being spatially structured is fulfilled, the test
#' is performed and the function returns an object of class \code{randtest} containing the
#' results of the test. The observed value (\code{obs}) is the adjusted R-squared of the
#' SSEF. Otherwise, a message is printed and \code{NULL} is returned invisibly.
#'
#' @author David Bauman and Jason Vleminckx, \email{davbauman@@gmail.com},
#' \email{jasvlx86@@gmail.com}
#'
#' @seealso \code{\link[vegan]{varpart}}, \code{\link{listw.select}}, \code{\link{listw.candidates}}, \code{\link{mem.select}}
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
#' Borcard D., Legendre P. and Drapeau P. (1992) Partialling out the spatial component of
#' ecological variation. Ecology, 73(3), 1045--1055
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
#' \donttest{
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
#' modsel.Y <- listw.select(Y, candidates, method = "FWD", MEM.autocor = "positive",
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
#' # Note: since vegan 2.6-4, the shared fraction (SSEF) is labelled [c] in the output of
#' # varpart(), and is the third row of VP$part$indfract.
#'
#' # Test of the shared space-environment fraction (SSEF):
#' SSEF.test <- envspace.test(Y, env, coord, MEM.spe,
#'                            listw.env = candidates[[modsel.env$best.id]],
#'                            regular = FALSE, nperm = 999)
#' SSEF.test
#'
#' # The SSEF is highly significant, indicating a potential induced spatial dependence.
#' }
#' }
#' @importFrom vegan rda anova.cca RsquareAdj
#' @importFrom stats runif
#' @export

"envspace.test" <- function(spe,
                            env,
                            coord,
                            MEM.spe,
                            listw.env,
                            MEM.autocor = c("positive", "negative", "all"),
                            regular = FALSE,
                            nperm = 999,
                            MSR.method = "singleton",
                            alpha = 0.05) {

  # 1. Checks and formatting of the inputs ----
  # **************************************************************************
  # **************************************************************************

  MEM.autocor <- match.arg(MEM.autocor)

  # All data tables are handled as numeric matrices from here on (a vector becomes a
  # one-column matrix), so that the same code applies to univariate and multivariate cases
  spe     <- as.matrix(spe)
  env     <- as.matrix(env)
  MEM.spe <- as.matrix(MEM.spe)
  coord   <- as.matrix(coord)
  n       <- nrow(coord)

  if (any(is.na(spe)) | any(is.na(coord)) | any(is.na(env)) | any(is.na(MEM.spe)))
    stop("NA entries in spe, env, coord or MEM.spe")
  if (!is.numeric(env))
    stop("'env' must only contain numeric variables")
  if (nrow(spe) != n | nrow(env) != n | nrow(MEM.spe) != n)
    stop("spe, env, coord and MEM.spe must have the same number of rows")
  if (inherits(listw.env, what = "listw") == FALSE)
    stop("listw.env is not an object of class 'listw'")
  if (ncol(MEM.spe) >= n - 1)
    stop("n-1 MEM variables in 'MEM.spe'. Select a subset to avoid saturated model issues.")
  # The model including both env and MEM.spe must not be saturated either, otherwise its
  # adjusted R2 (and thus the SSEF) cannot be estimated
  if (ncol(env) + ncol(MEM.spe) >= n - 1)
    stop("ncol(env) + ncol(MEM.spe) >= n - 1: the model with both 'env' and 'MEM.spe' ",
         "is saturated. Reduce the number of explanatory variables.")

  # For the torus-translation test, the grid structure is derived (and checked) right
  # away, so that the function stops before the (longer) preliminary tests if the
  # sampling design is not regular
  if (isTRUE(regular)) grid <- .envspace.torus.grid(coord)

  # 2. Preliminary tests ----
  # **************************************************************************
  # **************************************************************************

  # 2.1. Does the environment display a significant spatial structure? ----
  # **************************************************************************
  global <- mem.select(env, listw.env, method = "global",
                       MEM.autocor = MEM.autocor, alpha = alpha)$global.test
  if (MEM.autocor != "all")
    pval <- global$pvalue
  else
    pval <- c(global$positive$pvalue, global$negative$pvalue)

  if (length(which(pval <= alpha)) == 0) {
    message("No significant spatial structure detected in 'env'; the test was not performed")
    return(invisible(NULL))
  }

  # 2.2. Does the response display a significant spatial structure? ----
  # **************************************************************************
  test.spe <- anova.cca(rda(spe, MEM.spe))
  if (test.spe$Pr[1] > alpha) {
    message("No significant spatial structure detected in 'spe'; the test was not performed")
    return(invisible(NULL))
  }

  # 2.3. Is the response significantly related to the environment? ----
  # **************************************************************************
  test.env <- anova.cca(rda(spe, env))
  if (test.env$Pr[1] > alpha) {
    message("No significant relation between 'spe' and 'env'; the test was not performed")
    return(invisible(NULL))
  }

  # 3. Observed SSEF ----
  # **************************************************************************
  # **************************************************************************

  # The SSEF is computed from its definition (fraction [b] of Borcard et al. 1992):
  #   SSEF = R2(env) + R2(MEM.spe) - R2(env + MEM.spe)
  # rather than extracted from the output of vegan::varpart() by position. The order and
  # labels of the fractions in varpart() changed in vegan 2.6-4 (the shared fraction moved
  # from the 2nd to the 3rd row of $part$indfract), which silently made earlier versions
  # of this function test the pure spatial fraction instead of the SSEF.

  # R2 of the spatial model: spe and MEM.spe never change across permutations, so this
  # quantity is computed once and reused for every null replicate
  R2.mem <- .envspace.R2(spe, MEM.spe)

  # Observed SSEF, both unadjusted and adjusted (named vector)
  SSEF.obs <- .envspace.ssef(spe, env, MEM.spe, R2.mem = R2.mem)

  # Direction of the test, defined from the sign of the unadjusted SSEF (a negative SSEF
  # indicates that env and MEM.spe have opposite effects on spe; see Bauman et al. 2018c)
  alternative <- ifelse(SSEF.obs[["r.squared"]] < 0, "less", "greater")

  # 4. Null distribution of the SSEF ----
  # **************************************************************************
  # **************************************************************************

  # 4.1. Null environmental datasets ----
  # **************************************************************************
  # A list of nperm null versions of env, each preserving the spatial structure of the
  # original environmental variables but breaking their relationship with spe

  if (isTRUE(regular)) {

    # 4.1.1. Torus translations (regular designs only) ----
    # ==========================================================================
    n.config <- 2 * grid$nx * grid$ny   # translations x (no reflection, point reflection)
    if (nperm >= n.config)
      message(sprintf(paste0("Only %d distinct torus translations (including point ",
                             "reflections) exist for this %d x %d design, so that many of ",
                             "the %d null replicates are duplicates (p-value resolution ",
                             "limited to about 1/%d)"),
                      n.config, grid$nx, grid$ny, nperm, n.config))
    env.null <- lapply(seq_len(nperm), function(k)
      env[.envspace.torus.draw(grid), , drop = FALSE])

  } else {

    # 4.1.2. Moran Spectral Randomisation (any sampling design) ----
    # ==========================================================================
    env.null <- msr(env, listw.env, nrepet = nperm, method = MSR.method,
                    simplify = FALSE)
  }

  # 4.2. SSEF of each null replicate ----
  # **************************************************************************
  SSEF.sim <- vapply(env.null, function(E)
    .envspace.ssef(spe, E, MEM.spe, R2.mem = R2.mem)[["adj.r.squared"]],
    numeric(1))

  # 5. Output ----
  # **************************************************************************
  # **************************************************************************

  return(as.randtest(sim = SSEF.sim, obs = SSEF.obs[["adj.r.squared"]],
                     alter = alternative, call = match.call()))
}


# Internal helpers (not exported) ----
# ****************************************************************************
# ****************************************************************************

# Unadjusted and adjusted R2 of the RDA of Y on X ----
# ****************************************************************************
# Returns a named numeric vector c(r.squared, adj.r.squared). The values are taken by
# name from vegan::RsquareAdj(), so they do not depend on any ordering convention.
.envspace.R2 <- function(Y, X) {
  unlist(RsquareAdj(rda(Y, X))[c("r.squared", "adj.r.squared")])
}

# Shared fraction of a two-table variation partitioning, from its definition ----
# ****************************************************************************
# SSEF = R2(X1) + R2(X2) - R2(X1 + X2), here with X1 = env and X2 = MEM, computed both
# for unadjusted and adjusted R2. 'R2.mem' (the R2 of Y on MEM) can be supplied to avoid
# recomputing it at each permutation.
.envspace.ssef <- function(Y, env, MEM, R2.mem = .envspace.R2(Y, MEM)) {
  .envspace.R2(Y, env) + R2.mem - .envspace.R2(Y, cbind(env, MEM))
}

# Grid structure of a regular sampling design, for torus translations ----
# ****************************************************************************
# Converts the coordinates into integer grid positions (column 'ix', row 'iy'), and stops
# if the design is not a complete grid (or transect) of evenly spaced sites. Grid
# positions are derived from the coordinate values, so the order of the rows and the units
# of the coordinates do not matter.
.envspace.torus.grid <- function(coord) {

  coord <- as.matrix(coord)
  if (ncol(coord) == 1) coord <- cbind(coord, 0)   # transect given as a single column
  if (ncol(coord) != 2)
    stop("'coord' must have one (transect) or two (grid) columns for the ",
         "torus-translation test (regular = TRUE)")

  # Sorted unique coordinate values along each axis
  ux <- sort(unique(coord[, 1]))
  uy <- sort(unique(coord[, 2]))

  # Evenly spaced axis: all the gaps between consecutive coordinate values are equal
  # (up to a numerical tolerance relative to the extent of the axis)
  is.even <- function(u) {
    if (length(u) < 3) return(TRUE)
    gaps <- diff(u)
    diff(range(gaps)) <= 1e-8 * max(1, diff(range(u)))
  }
  if (!is.even(ux) | !is.even(uy))
    stop("The sites are not evenly spaced along both axes: the sampling design is not ",
         "regular and the torus-translation test cannot be used. Use regular = FALSE (MSR).")

  # Integer grid positions of each site
  ix <- match(coord[, 1], ux)
  iy <- match(coord[, 2], uy)
  nx <- length(ux)
  ny <- length(uy)

  # Complete grid: every grid node is sampled exactly once
  if (nx * ny != nrow(coord) | anyDuplicated(cbind(ix, iy)) > 0)
    stop("The sampling design is not a complete regular grid (each grid node must be ",
         "sampled exactly once): the torus-translation test cannot be used. ",
         "Use regular = FALSE (MSR).")

  # Look-up matrix: index of the site (row of the data) located at each grid node
  cell <- matrix(NA_integer_, nrow = nx, ncol = ny)
  cell[cbind(ix, iy)] <- seq_len(nrow(coord))

  list(ix = ix, iy = iy, nx = nx, ny = ny, cell = cell)
}

# One random torus translation ----
# ****************************************************************************
# Returns a permutation 'idx' of the rows such that env[idx, ] is a torus translation of
# env: the whole environmental map is shifted by a random number of cells along each axis
# (uniform over the whole extent of the axis, with wrap-around), and is rotated by 180
# degrees (point reflection) with probability 0.5. Both operations preserve the spatial
# autocorrelation structure of env.
.envspace.torus.draw <- function(grid) {

  # Random shifts, in number of grid cells (0 to nx - 1, and 0 to ny - 1)
  sx <- sample.int(grid$nx, 1) - 1
  sy <- sample.int(grid$ny, 1) - 1

  # New grid position of the environmental values of each site after the shift
  new.ix <- (grid$ix - 1 + sx) %% grid$nx + 1
  new.iy <- (grid$iy - 1 + sy) %% grid$ny + 1

  # Point reflection (180-degree rotation) of the whole grid, with probability 0.5
  if (runif(1) < 0.5) {
    new.ix <- grid$nx + 1 - new.ix
    new.iy <- grid$ny + 1 - new.iy
  }

  # dest[i] is the site that receives the environmental values of site i, so that the
  # environmental values of site j in the null dataset are those of site order(dest)[j]
  dest <- grid$cell[cbind(new.ix, new.iy)]
  order(dest)
}
