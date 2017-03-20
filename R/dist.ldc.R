#' Dissimilarity matrices for community composition data
#' 
#' Compute dissimilarity indices for ecological data matrices. The dissimilarity
#' indices computed by this function are those described in Legendre & De 
#' Cáceres (2013). In the name of the function, 'ldc' stands for the author's 
#' names. Twelve of these 21 indices are not readily available in other R 
#' package functions; four of them can, however, be computed in two computation 
#' steps in \code{vegan}.
#' 
#' @param Y Community composition data. The object class can be either 
#'   \code{data.frame} or \code{matrix}.
#'   
#' @param method One of the 21 dissimilarity coefficients available in the 
#'   function: \code{"hellinger"}, \code{"chord"}, \code{"chisquare"}, 
#'   \code{"profiles"}, \code{"percentdiff"}, \code{"ruzicka"}, 
#'   \code{"divergence"}, \code{"canberra"}, \code{"whittaker"}, 
#'   \code{"wishart"}, \code{"kulczynski"}, \code{"jaccard"}, \code{"sorensen"},
#'   \code{"ochiai"}, \code{"ab.jaccard"}, \code{"ab.sorensen"}, 
#'   \code{"ab.ochiai"}, \code{"ab.simpson"}, \code{"euclidean"}, 
#'   \code{"manhattan"}, \code{"modmeanchardiff"}. See Details. Names can be 
#'   abbreviated to a non-ambiguous set of first letters. Default: 
#'   \code{method="hellinger"}.
#'   
#' @param binary If \code{binary=TRUE}, the data are transformed to 
#'   presence-absence form before computation of the dissimilarities. Default 
#'   value: \code{binary=FALSE}, except for the Jaccard, Sørensen and Ochiai 
#'   indices where \code{binary=TRUE}.
#'   
#' @param samp If \code{samp=TRUE}, the abundance-based distances (ab.jaccard, 
#'   ab.sorensen, ab.ochiai, ab.simpson) are computed for sample data. If 
#'   \code{samp=FALSE}, binary indices are computed for true population data.
#'   
#' @param silent If \code{silent=FALSE}, informative messages sent to users will
#'   be printed to the R console. Use \code{silent=TRUE} is called on a 
#'   numerical simulation loop, for example.
#'   
#' @details The dissimilarities computed by this function are the following. 
#'   Indices i and k designate two rows (sites) of matrix Y, j designates a 
#'   column (species). D[ik] is the dissimilarity between rows i and k. p is the
#'   number of columns (species) in Y; pp is the number of species present in 
#'   one or the other site, or in both. y[i+] is the sum of values in row i; 
#'   same for y[k+]. y[+j] is the sum of values in column j. y[++] is the total 
#'   sum of values in Y. The indices are computed by functions written in C for 
#'   greater computation speed with large data matrices. \itemize{ \item Group 1
#'   - D computed by transformation of Y followed by Euclidean distance 
#'   \itemize{\item Hellinger D, D[ik] = 
#'   sqrt(sum((sqrt(y[ij]/y[i+])-sqrt(y[kj]/y[k+]))^2)) \item chord D, D[ik] = 
#'   sqrt(sum((y[ij]/sqrt(sum(y[ij]^2))-y[kj]/sqrt(sum(y[kj]^2)))^2)) \item 
#'   chi-square D, D[ik] = sqrt(y[++] sum((1/j[+j])(y[ij]/y[i+]-y[kj]/y[k+])^2))
#'   \item species profiles D, D[ik] = sqrt(sum((y[ij]/y[i+]-y[kj]/y[k+])^2)) }
#'   
#'   \item Group 2 - Other D functions appropriate for beta diversity studies 
#'   where A = sum(min(y[ij],y[kj])), B = y[i+]-A, C = y[k+]-A \itemize{\item 
#'   percentage difference D, D[ik] = (sum(abs(y[ij]-y[k,j])))/(y[i+]+y[k+]) or 
#'   else, D[ik] = (B+C)/(2A+B+C) = \item Ružička D, D[ik] = 
#'   1-(sum(min(y[ij],y[kj])/sum(max(y[ij],y[kj])) or else, D[ik] = 
#'   (B+C)/(A+B+C) \item coeff. of divergence D, D[ik] = 
#'   sqrt((1/pp)sum(((y[ij]-y(kj])/(y[ij]+y(kj]))^2)) \item Canberra metric D, 
#'   D[ik] = (1/pp)sum(abs(y[ij]-y(kj])/(y[ij]+y(kj])) \item Whittaker D, D[ik] 
#'   = 0.5*sum(abs(y[ij]/y[i+]-y(kj]/y[k+])) \item Wishart D, D[ik] = 
#'   1-sum(y[ij]y[kj])/(sum(y[ij]^2)+sum(y[kj]^2)-sum(y[ij]y[kj])) \item 
#'   Kulczynski D, D[ik] = 
#'   1-0.5((sum(min(y[ij],y[kj])/y[i+]+sum(min(y[ij],y[kj])/y[k+])) }
#'   
#'   \item Group 3 - Classical indices for binary data; they are appropriate for
#'   beta diversity studies. Value a is the number of species found in both i 
#'   and k, b is the number of species in site i not found in k, and c is the 
#'   number of species found in site k but not in i. The D matrices are 
#'   square-root transformed, as in dist.binary of ade4; the user-oriented 
#'   reason for this transformation is explained below. \itemize{\item Jaccard 
#'   D, D[ik] = sqrt((b+c)/(a+b+c)) \item Sørensen D, D[ik] = 
#'   sqrt((b+c)/(2a+b+c)) \item Ochiai D, D[ik] = sqrt(1 - a/sqrt((a+b)(a+c))) }
#'   
#'   \item Group 4 - Abundance-based indices of Chao et al. (2006) for 
#'   quantitative abundance data. These functions correct the index for species 
#'   that have not been observed due to sampling errors. For the meaning of the 
#'   U and V notations, see Chao et al. (2006, section 3). When 
#'   \code{samp=TRUE}, the abundance-based distances (ab.jaccard, ab.sorensen, 
#'   ab.ochiai, ab.simpson) are computed for sample data. If \code{samp=FALSE}, 
#'   indices are computed for true population data. – Do not use indices of 
#'   group 4 with \code{samp=TRUE} on presence-absence data; the indices are not
#'   meant to accommodate this type of data. If \code{samp=FALSE} is used with 
#'   presence-absence data, the indices are the regular {Jaccard, Sorensen, 
#'   Ochiai, Simpson} indices. On output, however, the D matrices are not 
#'   square-rooted, contrary to the {Jaccard, Sorensen, Ochiai} indices in 
#'   section 3 which are square-rooted. \itemize{\item abundance-based Jaccard 
#'   D, D[ik] = 1-(UV/(U+V-UV)) \item abundance-based Sørensen D, D[ik] = 
#'   1-(2UV/(U+V)) \item abundance-based Ochiai D, D[ik] = 1-sqrt(UV) \item 
#'   abundance-based Simpson D, D[ik] = 1-(UV/(UV+min((U-UV),(V-UV)))) }
#'   
#'   \item Group 5 – General-purpose dissimilarities that do not have an upper 
#'   bound (maximum D value). They are inappropriate for beta diversity studies.
#'   \itemize{\item Euclidean D, D[ik] = sqrt(sum(y[ij]-y[kj])^2) \item 
#'   Manhattan D, D[ik] = sum(abs(y[ij] - y[ik])) \item modified mean character 
#'   difference, D[ik] = (1/pp) sum(abs(y[ij] - y[ik])) } } The properties of 
#'   all dissimilarities available in this function (except Ružička D) were 
#'   described and compared in Legendre & De Cáceres (2013), who showed that 
#'   most of these dissimilarities are appropriate for beta diversity studies. 
#'   Inappropriate are the Euclidean, Manhattan, modified mean character 
#'   difference, species profile and chi-square distances. Most of these 
#'   dissimilarities have a maximum value of either 1 or sqrt(2). Three 
#'   dissimilarities (Euclidean, Manhattan, Modified mean character difference) 
#'   do not have an upper bound and are thus inappropriate for beta diversity 
#'   studies. The chi-square distance has an upper bound of
#'   sqrt(2*(sum(Y))).\cr\cr The Euclidean, Hellinger, chord, chi-square and
#'   species profiles dissimilarities have the property of being Euclidean,
#'   meaning that they never produce negative eigenvalues in principal
#'   coordinate analysis. The Canberra, Whittaker, percentage difference,
#'   Wishart and Manhattan coefficients are Euclidean when they are square-root
#'   transformed (Legendre & De Cáceres 2013, Table 2). The distance forms (1-S)
#'   of the Jaccard, Sørensen and Ochiai similarity (S) coefficients are
#'   Euclidean after taking the square root of (1-S) (Legendre & Legendre 2012,
#'   Table 7.2). The D matrices resulting from these three coefficients are
#'   outputted in the form sqrt(1-S), as in function \code{dist.binary} of ade4,
#'   because that form is Euclidean and will thus produce no negative
#'   eigenvalues in principal coordinate analysis. \cr\cr The Hellinger, chord,
#'   chi-square and species profile dissimilarities are computed using the
#'   two-step procedure developed by Legendre & Gallagher (2001). The data are
#'   first transformed using either the row marginals, or the row and column
#'   marginals in the case of the chi-square distance. The dissimilarities are
#'   then computed from the transformed data using the Euclidean distance
#'   formula. As a consequence, these four dissimilarities are necessarily
#'   Euclidean. D matrices for other binary coefficients can be computed in two
#'   ways: either by using function \code{dist.binary} of ade4, or by choosing
#'   option \code{binary=TRUE}, which transforms the abundance data to binary
#'   form, and using one of the quantitative indices of the present function.
#'   Table 1 of Legendre & De Cáceres (2013) shows the incidence-based
#'   (presence-absence-based) indices computed by the various indices using
#'   binary data. \cr\cr The Euclidean distance computed on untransformed
#'   presence-absence or abundance data produces non-informative and incorrect
#'   ordinations, as shown in Legendre & Legendre (2012, p. 300) and in Legendre
#'   & De Cáceres (2013). However, the Euclidean distance computed on
#'   log-transformed abundance data produces meaningful ordinations in principal
#'   coordinate analysis (PCoA). Nonetheless, it is easier to compute a PCA of
#'   log-transformed abundance data instead of a PCoA; the resulting ordination
#'   with scaling 1 will be meaningful. Messages are printed to the R console
#'   indicating the Euclidean status of the computed dissimilarity matrices.
#'   Note that for the chi-square distance, the columns that sum to zero are
#'   eliminated before calculation of the distances, thus preventing divisions
#'   by zero in the calculation of the chi-square transformation.
#'   
#' @return A dissimilarity matrix, with class \code{dist}.
#'   
#' @references Chao, A., R. L. Chazdon, R. K. Colwell and T. J. Shen. 2006. 
#'   Abundance-based similarity indices and their estimation when there are 
#'   unseen species in samples. Biometrics 62: 361–371.
#'   
#'   Legendre, P. and M. De Cáceres. 2013. Beta diversity as the variance of 
#'   community data: dissimilarity coefficients and partitioning. Ecology 
#'   Letters 16: 951-963.
#'   
#'   Legendre, P. and E. D. Gallagher, E.D. 2001. Ecologically meaningful 
#'   transformations for ordination of species data. Oecologia 129: 271–280.
#'   
#'   Legendre, P. and Legendre, L. 2012. Numerical Ecology. 3rd English edition.
#'   Elsevier Science BV, Amsterdam.
#'   
#' @author Pierre Legendre \email{pierre.legendre@@umontreal.ca} and Naima Madi
#'   
#' @examples 
#' 
#' if(require("vegan", quietly = TRUE)) {
#' data(mite)
#' mat1  = as.matrix(mite[1:10, 1:15])   # No column has a sum of 0
#' mat2 = as.matrix(mite[61:70, 1:15])   # 7 of the 15 columns have a sum of 0
#' 
#' #Example 1: compute Hellinger distance for mat1
#' D.out = dist.ldc(mat1,"hellinger")
#' 
#' #Example 2: compute chi-square distance for mat2
#' D.out = dist.ldc(mat2,"chisquare")
#' 
#' #Example 3: compute percentage difference dissimilarity for mat2
#' D.out = dist.ldc(mat2,"percentdiff")
#' 
#' }
#' 
#' 
#' @export dist.ldc
dist.ldc <-
    function(Y,
        method = "hellinger",
        binary = FALSE,
        samp = TRUE,
        silent = FALSE)
    {
        epsilon <- sqrt(.Machine$double.eps)
        method <-
            match.arg(
                method,
                c(
                    "hellinger",
                    "chord",
                    "chisquare",
                    "profiles",
                    "percentdiff",
                    "ruzicka",
                    "divergence",
                    "canberra",
                    "whittaker",
                    "wishart",
                    "kulczynski",
                    "jaccard",
                    "sorensen",
                    "ochiai",
                    "ab.jaccard",
                    "ab.sorensen",
                    "ab.ochiai",
                    "ab.simpson",
                    "euclidean",
                    "manhattan",
                    "modmeanchardiff"
                )
            )
        #
        Y <- as.matrix(Y)
        n <- nrow(Y)
        
        if ((n == 2) &
                (dist(Y)[1] < epsilon))
            stop("Y only contains two rows and they are identical")
        if (any(Y < 0))
            stop("Data contain negative values\n")
        if (any(is.na(Y)))
            stop("Data contain 'NA' values\n")
        # Transform data to presence-absence before computing the binary coefficients
        # or if the user choses binary=TRUE
        if (any(method == c("jaccard", "sorensen", "ochiai")))
            binary = TRUE
        if (binary)
            Y <- ifelse(Y > 0, 1, 0)
        #
        # Group 1 indices
        switch(
            method,
            hellinger = {
                YY = .Call("transform_mat", Y, "hellinger")
                D = .Call("euclidean", YY)
                if (!silent)
                    cat("Info -- This coefficient is Euclidean***\n")
            },
            chord = {
                YY = .Call("transform_mat", Y, "chord")
                D = .Call("euclidean", YY)
                if (!silent)
                    cat("Info -- This coefficient is Euclidean\n")
            },
            chisquare = {
                YY = .Call("transform_mat", Y, "chisquare")
                D = .Call("euclidean", YY)
                if (!silent)
                    cat("Info -- This coefficient is Euclidean\n")
            },
            profiles = {
                YY = .Call("transform_mat", Y, "profiles")
                D = .Call("euclidean", YY)
                if (!silent)
                    cat("Info -- This coefficient is Euclidean\n")
                
                # Group 2 indices
            },
            percentdiff = {
                D = .Call("percentdiff", Y)
                if (!silent)
                    cat("Info -- For this coefficient, sqrt(D) would be Euclidean\n")
            },
            ruzicka = {
                D = .Call("ruzicka", Y)
                if (!silent)
                    cat("Info -- For this coefficient, sqrt(D) would be Euclidean\n")
            },
            divergence = {
                D = .Call("divergence", Y)
                if (!silent)
                    cat("Info -- This coefficient is Euclidean\n")
            },
            canberra = {
                D = .Call("canberra", Y)
                if (!silent)
                    cat("Info -- For this coefficient, sqrt(D) would be Euclidean\n")
            },
            whittaker = {
                D = .Call("whittaker", Y)
                if (!silent)
                    cat("Info -- For this coefficient, sqrt(D) would be Euclidean\n")
            },
            wishart = {
                D = .Call("wishart", Y)
                if (!silent)
                    cat("Info -- For this coefficient, sqrt(D) would be Euclidean\n")
            },
            kulczynski = {
                D = .Call("kulczynski", Y)
                if (!silent)
                    cat("Info -- This coefficient is not Euclidean\n")
                
                # Group 3 indices
            },
            jaccard = {
                D = .Call("binary_D", Y, "jaccard")
                if (!silent)
                    cat("Info -- D is Euclidean because the function outputs D[jk] = sqrt(1-S[jk])\n")
            },
            sorensen = {
                D = .Call("binary_D", Y, "sorensen")
                if (!silent)
                    cat("Info -- D is Euclidean because the function outputs D[jk] = sqrt(1-S[jk])\n")
            },
            ochiai = {
                D = .Call("binary_D", Y, "ochiai")
                if (!silent)
                    cat("Info -- D is Euclidean because the function outputs D[jk] = sqrt(1-S[jk])\n")
                
                # Group 4 indices
            },
            ab.jaccard = {
                D = .Call("chao_C", Y, "Jaccard", samp)
                if (!silent) {
                    cat("Info -- This coefficient is not Euclidean\n")
                    if (!samp) {
                        cat("If data are presence-absence, sqrt(D) will be Euclidean\n")
                    }
                }
            },
            ab.sorensen = {
                D = .Call("chao_C", Y, "Sorensen", samp)
                if (!silent) {
                    cat("Info -- This coefficient is not Euclidean\n")
                    if (!samp) {
                        cat("If data are presence-absence, sqrt(D) will be Euclidean\n")
                    }
                }
            },
            ab.ochiai = {
                D = .Call("chao_C", Y, "Ochiai", samp)
                if (!silent) {
                    cat("Info -- This coefficient is not Euclidean\n")
                    if (!samp) {
                        cat("If data are presence-absence, sqrt(D) will be Euclidean\n")
                    }
                }
            },
            ab.simpson = {
                D = .Call("chao_C", Y, "Simpson", samp)
                if (!silent) {
                    cat("Info -- This coefficient is not Euclidean\n")
                    if (!samp) {
                        cat("If data are presence-absence, sqrt(D) will be Euclidean\n")
                    }
                }
                
                # Group 5 indices
            },
            euclidean = {
                D = .Call("euclidean", Y)
                if (!silent) {
                    cat("Info -- This coefficient is Euclidean\n")
                    cat("Info -- This coefficient does not have an upper bound (no fixed D.max)\n")
                }
            },
            manhattan = {
                D = .Call("manhattan", Y)
                if (!silent) {
                    cat("Info -- For this coefficient, sqrt(D) would be Euclidean\n")
                    cat("Info -- This coefficient does not have an upper bound (no fixed D.max)\n")
                }
            },
            modmeanchardiff = {
                D = .Call("modmean", Y)
                if (!silent) {
                    cat("Info -- This coefficient is not Euclidean\n")
                    cat("Info -- This coefficient does not have an upper bound (no fixed D.max)\n")
                }
            }
        )
        #
        D <- as.dist(D)
    }
