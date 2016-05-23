#' Compute LCBD from any D matrix
#'
#' Compute LCBD indices (Legendre & De Cáceres 2013) from a symmetric
#' dissimilarity matrix (D) or from a beta component matrix (Repl, RichDiff or
#' AbDiff, or Nes) (Legendre 2014).
#'
#' @param D A dissimilarity or beta diversity component matrix (class
#'   \code{dist} or \code{matrix}.
#' @param sqrt.D Take the square root of the dissimilarities in matrix D before
#'   computing the LCBD indices.
#' @param save.D If \code{save.D} is \code{TRUE}, the dissimilarity matrix will
#'   appear in the output list.
#'
#' @details Use \code{sqrt.D = TRUE} when computing LCBD indices for most of the
#'   replacement and richness/abundance difference indices computed by function
#'   \code{beta.div.comp}, as well as for the corresponding D matrices. See
#'   Table S1.4 in Appendix S1 of Legendre (2014) to identify the matrices that
#'   are Euclidean without taking the square root of the individual values. Only
#'   the RichDiffS (for presence-absence data) and AbDiff%diff indices (for
#'   abundance data) of the Sørensen group in the Podani family have that
#'   property. In all other cases, use \code{sqrt.D = TRUE}. When computing LCBD
#'   from a D matrix, use \code{sqrt = TRUE} if the D matrix is not Euclidean.
#'   The Euclidean property can be checked with function \code{is.euclid} of
#'   \code{ade4}. BDtotal statistics are comparable among data sets having the
#'   same or different numbers of sampling units (n), provided that the sampling
#'   units are of the same size or represent the same sampling effort and that
#'   BDtotal is computed with the same D index. Function \code{LCBD.comp}
#'   produces the same (SStotal, BDtotal, LCBD) results as function
#'   \code{beta.div}. Note, however, that the latter produces other interesting
#'   results (p.LCBD, SCBD). Function \code{LCBD.comp} should then only be used
#'   to compute LCBD indices from dissimilarity matrices that cannot be computed
#'   by function \code{beta.div}, e.g. genetic D matrices, or from replacement
#'   and richness difference matrices produced by function \code{beta.div.comp}.
#'   Significance of the LCBD indices cannot be tested when their calculation
#'   starts from a D matrix because the testing procedure involves permutation
#'   of the columns of raw data.
#'
#' @return A list containing the following results:
#' \itemize{ \item \code{SStotal_BDtotal}: Total sum of squares and total beta diversity [= Var(Y)] of
#'   the data matrix. \item \code{LCBD}: Vector of Local contributions to beta diversity
#'   (LCBD) for the sites. \item \code{D}: The input dissimilarity matrix, class \code{dist}; only
#'   if \code{save.D=TRUE}}.
#'
#' @author Pierre Legendre \email{pierre.legendre@@umontreal.ca}
#'
#' @references Legendre, P. 2014. Interpreting the replacement and richness
#'   difference components of beta diversity. Global Ecology and Biogeography
#'   23: 1324-1334.
#'
#'   Legendre, P. & M. De Cáceres. 2013. Beta diversity as the variance of
#'   community data: dissimilarity coefficients and partitioning. Ecology
#'   Letters 16: 951-963.
#'
#' @keywords spatial
#' @examples
#'
#' ### Example 1
#' ### Compute the Hellinger distance, then the LCBD indices.
#' if(require("vegan", quietly = TRUE)){
#' data(mite)
#' mite.hel = decostand(mite, "hellinger")
#' mite.D = dist(mite.hel)
#' out.mite.D = LCBD.comp(mite.D, sqrt.D=FALSE)
#' }
#'
#' ### Example 2
#' if(require("ade4", quietly = TRUE) & require("adegraphics", quietly = TRUE)){
#' data(doubs)
#' fish.sp = doubs$fish[-8,]   # Fish data; site 8 is removed because no fish were caught
#'
#' out.comp = beta.div.comp(fish.sp, coef="S", quant=TRUE)
#'
#' out.fish.D = LCBD.comp(out.comp$D, sqrt.D=TRUE)   # out.comp.D is not Euclidean
#' out.fish.D$SStotal_BDtotal
#' out.fish.Repl = LCBD.comp(out.comp$repl, sqrt.D=TRUE)   # out.comp$repl is not Euclidean
#' out.fish.Repl$SStotal_BDtotal
#' out.fish.AbDiff = LCBD.comp(out.comp$rich, sqrt.D=FALSE)   # out.comp$rich is Euclidean
#' out.fish.AbDiff$SStotal_BDtotal
#'
#' ### Plot maps of the LCBD indices
#' fish.xy = doubs$xy[-8,]   # Geographic coordinates; site 8 removed because no fish were caught
#'
#' # Map of LCBD indices for %difference dissimilarity
#' s.value(fish.xy, out.fish.D$LCBD, method="size", symbol = "circle",
#' col = c("white", "brown"), main = "Doubs fish LCBD, %difference D")
#'
#' # Map of LCBD indices for replacement component of %difference dissimilarity
#' s.value(fish.xy, out.fish.Repl$LCBD, method="size", symbol = "circle",
#' col = c("white", "brown"), main = "Doubs fish replacement LCBD")
#'
#' # Map of LCBD indices for abundance difference component of %difference dissimilarity
#' s.value(fish.xy, out.fish.AbDiff$LCBD, method="size", symbol = "circle", 
#' col = c("white", "brown"), main = "Doubs fish abundance diff. LCBD")
#' }
#'
#' @importFrom stats as.dist
#' @export LCBD.comp
#'

LCBD.comp <- function(D, sqrt.D = TRUE, save.D = FALSE) {
    #
    # Description --
    #
    # Compute LCBD indices (Legendre and De Cáceres 2013) from a dissimilarity
    # matrix (D) or beta div. component matrices (Repl, RichDiff or AbDiff, or Nes).
    #
    # Arguments --
    #
    # D : Dissimilarity or beta diversity component matrix, class=dist.
    # sqrt.D : Take sqrt() of components before computing LCBD.comp. Use
    #     sqrt.D=TRUE for the replacement and richness/abundance difference indices
    #     computed by beta.div.comp(), as well as for the corresponding D matrices.
    # =>  When computing LCBD from a D matrix, use sqrt=TRUE if the D matrix is not
    #     Euclidean. That property can be checked with function is.euclid() of ade4.
    #
    # Reference --
    #
    # Legendre, P. & De Cáceres, M. (2013) Beta diversity as the variance of
    # community data: dissimilarity coefficients and partitioning. Ecology
    # Letters 16: 951–963.
    #
    # License: GPL-2
    # Author:: Pierre Legendre, August 2013
    
    ### Internal function   # Legendre & Legendre 2012, eq. 9.42
    centre <- function(D, n)
        # Centre a square matrix D by matrix algebra
        # mat.cen = (I - 11'/n) D (I - 11'/n)
    {
        One <- matrix(1, n, n)
        mat <- diag(n) - One / n
        mat.cen <- mat %*% D %*% mat
    }
    ###
    D <- as.dist(D)
    x <- as.matrix(D)
    n <- nrow(x)
    
    if (sqrt.D) {
        # D is used in Gower centring
        SStotal <- sum(D) / n
        BDtotal <- SStotal / (n - 1)
        G <- centre(as.matrix(-0.5 * x), n)
    } else {
        # D^2 is used in Gower centring
        SStotal <- sum(D ^ 2) / n      # eq. 8
        BDtotal <- SStotal / (n - 1)   # eq. 3
        G <- centre(as.matrix(-0.5 * x ^ 2), n)
    }
    LCBD <-
        diag(G) / SStotal   # Legendre & De Caceres (2013), eq. 10b
    if (save.D) {
        out <- list(
            SStotal_BDtotal = c(SStotal, BDtotal),
            LCBD = LCBD,
            D = D
        )
    } else {
        out <- list(SStotal_BDtotal = c(SStotal, BDtotal),
                    LCBD = LCBD)
    }
}