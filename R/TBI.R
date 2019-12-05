#' TBI: Difference between multivariate observations at T1 and T2
#'
#' The function computes and tests Temporal Beta-diversity Indices (TBI) between
#' multivariate observations (frequency or presence-absence data) forming pairs
#' observed at time 1 (T1) and time 2 (T2). The data matrices may contain
#' abundance or presence-absence data, or other types of frequency-like data
#' (e.g. biomasses). TBI are dissimilarity indices that measure beta
#' differentiation through time. The indices are computed between T1 and T2 for
#' each site. The difference between species (or abundances-per-species) gains
#' (C/den) and losses (B/den) can be printed out and tested for significance.
#'
#' @param mat1,mat2 Two multivariate community composition or gene frequency
#'   data matrices (class data.frame or matrix) with the same number of rows and
#'   columns. The rows must correspond to the same objects (e.g. sites) and the
#'   columns to the same variables, e.g. species or alleles.
#'
#' @param method One of the following dissimilarity coefficients:
#'   \code{"\%difference"} (aka Bray-Curtis), \code{"ruzicka"}, \code{"chord"},
#'   \code{"hellinger"}, \code{"log.chord"}, \code{"sorensen"},
#'   \code{"jaccard"}, \code{"ochiai"}, \code{"euclidean"}. See Details. Names
#'   can be abbreviated to a non-ambiguous set of first letters. Default:
#'   \code{method="\%difference"}.
#'
#' @param pa.tr If \code{pa.tr=TRUE}, the data are transformed to binary (i.e.
#'   presence-absence, or pa) form. If \code{pa.tr=FALSE}, they are not.
#'
#' @param nperm Number of permutations for the tests of significance of the
#'   temporal beta indices and the permutation test of the B-C difference. Use
#'   \code{nperm} = 999 or 9999 in real studies.
#'
#' @param BCD If \code{BCD=TRUE}, the B and C components of the
#'   percentage difference (\code{method="\%difference"}) and Ruzicka
#'   (\code{method"ruzicka"}) indices are computed and presented in an output
#'   table with three columns: B/den, C/den, D=(B+C)/den, where den is the
#'   denominator of the index, i.e. (2A+B+C) for the percentage difference index
#'   and (A+B+C) for the Ruzicka index. See Details and Value. If
#'   \code{pa.tr=TRUE}, the B and C components are the numbers of species lost
#'   or gained, and D is either the Sorensen or the Jaccard dissimilarity. In
#'   the BCD output table, column B contains B/den, C/den, D=(B+C)/den, as in
#'   the case of the percentage difference and Ruzicka indices.
#'   
#'   If \code{BCD=FALSE}, that table is not produced. No table can be computed for
#'   indices other than the Ruzicka and percentage difference or their binary
#'   forms. 
#'
#' @param replace If \code{FALSE} (default value),
#'   sampling is done without replacement, producing a regular permutation test.
#'   If \code{TRUE}, sampling is done with replacement for the
#'   test of significance; the method is then bootstrapping.
#'
#' @param test.BC If \code{TRUE}, the difference between
#'   species (or abundances-per-species) gains (C/den) and species (or
#'   abundances-per-species) losses (B/den) is tested in a parametric paired
#'   t-test computed by function \code{t.test} of \code{stats}. If
#'   \code{FALSE}, the test is not computed. 
#'
#' @param test.t.perm If \code{TRUE}, the difference
#'   between species (or abundances-per-species) gains (C/den) and species (or
#'   abundances-per-species) losses (B/den) is also tested in a permutational
#'   paired t-test computed by function \code{t.paired.perm}.If
#'   \code{FALSE}, the test is not computed.
#'
#' @param save.BC If \code{TRUE}, the original B and C
#'   values are returned in a matrix called \code{BC}, without division by den
#'   as in the output matrix \code{BCD.mat}. If \code{FALSE}, BC
#'   will get the value NA in the output list.
#'
#' @param seed. Seed for random number generator. If \code{NULL}, the random 
#'   number generator keeps going from the point it had reached in previous
#'   calculations. If \code{seed.} is an integer value instead of NULL,
#'   the random number generator is reset using that value. This allows users to
#'   repeat exactly a previous calculation launched with the same value of
#'   \code{seed.}; the sequence of generated random numbers will be exactly the
#'   same. 
#'
#' @param clock If \code{clock=TRUE}, the computation time is printed. This
#'   option is useful to predict the calculation time when \emph{n} and
#'   \emph{nperm} are large.
#'
#' @details For each object, the function tests the hypothesis (H0) that a
#' species assemblage is not exceptionally different between T1 and T2, compared
#' to assemblages that could have been observed at this site at T1 and T2 under
#' conditions corresponding to H0. If H0 is rejected, the object is recognized
#' as exceptionally different from the other objects for its difference between
#' T1 and T2.
#'
#' \emph{To fix ideas, an example in palaeoecology} - A researcher is studying
#' ancient and modern diatom communities in sediment cores. If a site displays
#' an exceptional difference between T1 and T2, the researcher is justified to
#' examine the reason for that difference. It could, for example, be caused by a
#' change in land use at that site, which has caused the difference to be larger
#' than at the other sites, compared to the differences caused by climate change
#' at all sites.
#'
#' The temporal beta diversity indices available in this function belong to four
#' groups, computed in different ways. \itemize{ \item Method
#' \code{"\%difference"} computes the percentage difference index, erroneously
#' called the Bray-Curtis index in some software packages ; it is the
#' quantitative form of the Sorensen index. Method \code{"ruzicka"} computes the
#' Ruzicka dissimilarity; this is one of the quantitative coefficients
#' corresponding to the Jaccard dissimilarity for binary data. When these
#' indices are used to compute ordinations by principal coordinate analysis, it
#' is recommended to take the square root of the dissimilarities before the
#' ordination analysis because these indices do not have the property of being
#' Euclidean. However, that precaution is not important here; the results of
#' permutation tests will be the same for these dissimilarities, square-rooted
#' or not. If \code{pa.tr=TRUE}, either the Sorensen or the Jaccard coefficient
#' are obtained by computing these two coefficients. \item Methods
#' \code{"chord"} (chord distance), \code{"hellinger"} (Hellinger distance) and
#' \code{"log.chord"} (log.chord distance) are obtained by transformation of the
#' species data, as described by Legendre & Borcard (2018), followed by
#' calculation of the Euclidean distance. For the Hellinger distance, the data
#' are square-rooted, then subjected to the chord transformation and the
#' Euclidean distance. For the log.chord distance, the data are transformed by
#' y' = log(y+1) using function log1p() of R, then subjected to the chord
#' transformation and the Euclidean distance. These three distances have the
#' Euclidean property (Legendre & Legendre 2012, Legendre & De Caceres 2013). If
#' \code{pa.tr=TRUE}, the  Ochiai distance for binary data,
#' sqrt(2)*sqrt(1-Ochiai similarity), is obtained from these three coefficients.
#' \item Methods \code{"jaccard"}, \code{"sorensen"}, \code{"ochiai"} implement
#' the Jaccard, Sorensen and Ochiai dissimilarities. For these coefficients, the
#' data are first transformed to presence-absence (\code{pa.tr} receives the
#' value \code{TRUE}), then the dissimilarities are computed using the
#' corresponding quantitative coefficients (Ruzicka, percentage difference, and
#' chord). \item The Euclidean distance is also available in this function. It
#' is not recommended for community composition or allele frequency data. One
#' can compute it for log-transformed abundance data that do not contain zeros,
#' or very few zeros (short gradients). }
#'
#' The temporal beta indices are tested for significance using permutation
#' tests. The hypotheses are the following: \itemize{ \item H0: the site under
#' study (e.g. a species assemblage) is not exceptionally different between T1
#' and T2, compared to assemblages that could have been observed at this site at
#' T1 and T2 under conditions corresponding to H0. The differences between T1
#' and T2 all belong to the same statistical population of differences. \item
#' H1: the site under study is exceptionally different between times T1 and T2.
#' }
#'
#' In the decomposition of the Ruzicka and percentage difference dissimilarities
#' or their presence-absence forms (Jaccard, Sorensen), the components B and C
#' are computed as follows: \itemize{ \item bj is the part of the abundance of
#' species j that is higher at time 1 than at time 2: bj = (y1j - y2j) if y1j >
#' y2j ; bj = 0 otherwise. B is the sum of the bj values for all species in the
#' group of species under study. It is the unscaled sum of species losses
#' between time 1 and time 2. In the BCD output table BCD.mat, column 1 contains
#' B/den where den is the denominator of the index, i.e. (2A+B+C) for the
#' percentage difference index and (A+B+C) for the Ruzicka index. \item cj is
#' the part of the abundance of species j that is higher at time 2 than at time
#' 1: cj = (y2j - y1j) if y2j > y1j ; cj = 0 otherwise. C is the sum of the cj
#' values for all species in the group of species under study. It is the
#' unscaled sum of species gains between time 1 and time 2. In the BCD output
#' table BCD.mat, column 2 contains C/den where den is the denominator of the
#' index, i.e. (2A+B+C) for the percentage difference index and (A+B+C) for the
#' Ruzicka index. }
#'
#' The original values of B and C for each site, without denominator, are also
#' available in the output table BC.
#'
#' \emph{Warning} - In real ecological studies, when the TBI test is applied to
#' data where some sites are highly impoverished due to pollution or other
#' extreme environmental situations, this situation may produce sites with very
#' few species (i.e. very low richness) and no species in common for the T1-T2
#' comparisons due to sampling variation at these impoverished sites. The TBI
#' dissimilarity will be high and the test may indicate a significant T1-T2
#' difference if most other sites have higher species richness. This would be a
#' correct statistical outcome for the test. When users of the method identify
#' sites showing significant TBI tests in data, they should check the species
#' richness of these sites at T1 and T2. Interpretation of the test results
#' should be done with caution when high and significant TBI indices are
#' associated with very low richness and no species in common between T1 and T2.
#'
#' @return Function TBI returns a list containing the following results:
#'   \itemize{ \item \code{TBI} The vector of Temporal Beta-diversity Indices
#'   (TBI) between observations at times T1 and T2 for each object.
#'
#'   \item \code{p.TBI} A corresponding vector of p-values. Significant p-values
#'   (e.g. p.TBI <= 0.05) indicate exceptional objects for the difference of
#'   their species composition.
#'
#'   \item \code{p.adj} The p-values are corrected for multiple testing using
#'   function p.adjust of {stats}. The adjustment is done using
#'   \code{method="holm"}, which is the default option of the \code{p.adjust}
#'   function.
#'
#'   \item \code{BCD.mat} An output table with four columns: B/den, C/den,
#'   D=(B+C)/den, and Change. The value den is the denominator of the index,
#'   i.e. (2A+B+C) for the percentage difference index and (A+B+C) for the
#'   Ruzicka index. The decomposition is such that D = B/den + C/den. Columns B
#'   and C indicate which of the D values are associated with large B (losses)
#'   or large C values (gains), before proceeding to the analysis and
#'   interpretation of the D values, using environmental or spatial explanatory
#'   variables, through regression or classification tree analysis. When B > C,
#'   the site has lost species or abundances-per-species between time 1 and time
#'   2; this is indicated by a "-" sign in column Change. On the contrary, if B
#'   < C, the site has gained species or abundances-per-species between time 1
#'   and time 2; this is indicated by a "+" sign in that column. Sites with
#'   equal amounts of losses and gains are marked with a "0". - The B/den and
#'   C/den values can be plotted in B-C plots, which are informative about the
#'   changes that occurred in the data set between the two surveys under study.
#'   - If \code{pa.tr} is TRUE, the B and C components are the numbers of
#'   spepcies losses and gains, and D is either the Sorensen or the Jaccard
#'   dissimilarity. - If \code{BCD=FALSE}, that table is not produced. No table
#'   is (or can be) computed for indices other than the Ruzicka and percentage
#'   difference indices or their binary forms.
#'
#'   \item \code{BCD.summary}	An output table with six columns: mean(B/den);
#'   mean(C/den); mean(D); B/(B+C) (which is mean(B/den) divided by mean(D));
#'   C/(B+C) (which is mean(C/den) divided by mean(D)). These values indicate,
#'   over all sites, which of the processes dominated (loss or gain of species
#'   or abundances-per-species) when site compositions changed between time 1
#'   and time 2. Change has the same meaning as in table \code{BCD.mat}; the
#'   sign indicates the direction of the mean change over all sites.
#'
#'   \item \code{t.test_B.C}	The results of a paired t-test (parametric) of
#'   significance of the difference between columns C/den and B/den of the
#'   \code{BCD.mat} table. If \code{test.t.perm=TRUE}, the difference between
#'   species gains (C/den) and losses (B/den) is also tested in a permutational
#'   paired t-test and the permutational p-value is shown in the output table.
#'   This result provides an overall test of the direction of change over all
#'   sites. It helps confirm the asymmetry between species (or
#'   abundances-per-species) gains (C/den) and species (or
#'   abundances-per-species) losses (B/den). A star in column p<=0.05 indicates
#'   a significant result of the parametric test at the 0.05 level.
#'
#'   \item \code{BC}	An output table with two columns: B and C. In this table,
#'   the B and C statistics are not divided by a denominator, contrary to the
#'   values B/den and C/den found in the output table \code{BCD.mat}. }
#'
#' @author Pierre Legendre \email{pierre.legendre@@umontreal.ca}
#'
#' @references Legendre, P. 2019. A temporal beta-diversity index to identify
#'   exceptional sites in space-time surveys. \emph{Ecology and Evolution} (in
#'   press).
#'
#'   Legendre, P. & M. De Caceres. 2013. Beta diversity as the variance of
#'   community data: dissimilarity coefficients and partitioning. \emph{Ecology
#'   Letters} 16: 951-963.
#'
#'   Legendre, P. & D. Borcard. 2018. Box-Cox-chord transformations for
#'   community composition data prior to beta diversity analysis.
#'   \emph{Ecography} 41: 1820-1824.
#'
#'   Legendre, P. & L. Legendre. 2012. \emph{Numerical Ecology. 3rd English
#'   edition.} Elsevier Science BV, Amsterdam.
#'
#'   van den Brink, P. J. & C. J. F. ter Braak. 1999. Principal response curves:
#'   analysis of time-dependent multivariate responses of biological community
#'   to stress. \emph{Environmental Toxicology and Chemistry} 18: 138-148.
#'
#' @seealso \code{\link{plot.TBI}}
#'
#' @examples
#' if(require("vegan", quietly = TRUE)) {
#'
#' ## Example 1 -
#'
#' ## Invertebrate communities subjected to insecticide treatment.
#'
#' ## As an example in their paper on Principal Response Curves (PRC method), van den
#' ## Brink & ter Braak (1999) used observations on the abundances of 178 invertebrate
#' ## species (macroinvertebrates and zooplankton) subjected to treatments in 12 mesocosms by
#' ## the insecticide chlorpyrifos. The mesocosms were sampled at 11 occasions. The data,
#' ## available in the {vegan} package, are log-transformed species abundances, ytranformed =
#' ## log(10*y+1).
#'
#' ## The data of survey #4 will be compared to those of survey #11 in this example.
#' ## Survey #4 was carried out one week after the insecticide treatment, whereas the fauna
#' ## of the mesocosms was considered by the authors to have fully recovered from the
#' ## insecticide treatment at survey #11.
#'
#' data(pyrifos)
#'
#' ## The mesocosms had originally been attributed at random to the treatments. However,
#' ## to facilitate presentation of the results, they will be listed here in order of
#' ## increased insecticide doses: {0, 0, 0, 0, 0.1, 0.1, 0.9, 0.9, 6, 6, 44, 44} micro g/L.
#'
#' ## Select the 12 data rows of surveys 4 and 11 from the data file and reorder them
#'
#' ord4 = c(38,39,41,47,37,44,40,46,43,48,42,45)
#'
#' ord11 = c(122,123,125,131,121,128,124,130,127,132,126,129)
#'
#' ## Run the TBI function
#'
#' res1 <- TBI(pyrifos[ord4,], pyrifos[ord11,], method = "%diff", nperm = 0, test.t.perm = FALSE)
#'
#' res1$BCD.mat
#'
#' ## Example 2 -
#'
#' ## This example uses the mite data available in vegan. Let us pretend that sites 1-20
#' ## represent T1 and sites 21-40 represent T2.
#'
#'
#' data(mite)
#'
#' # Run the TBI function
#'
#' res2 <- TBI(mite[1:20,], mite[21:40,], method = "%diff", nperm = 0, test.t.perm = FALSE)
#'
#' summary(res2)
#'
#' res2$BCD.mat
#'
#' }
#' @importFrom vegan decostand
#' @importFrom stats t.test
#' @export TBI
#'   

TBI <- function(mat1, mat2, method = "%difference", pa.tr = FALSE, nperm = 99, BCD = TRUE, replace = FALSE,
    test.BC = TRUE, test.t.perm = FALSE, save.BC = FALSE, seed. = NULL, clock = FALSE) {
    ### Internal functions --
    
    RuzickaD <- function(vec1, vec2, method = "ruzicka", BCD = FALSE, ref = TRUE) {
        #
        # Compute the Ruzicka dissimilarity (quantitative form of the Jaccard dissimilarity)
        # or the percentage difference (quantitative form of the Sorensen dissimilarity).
        # Algorithm applicable to two vectors only: data for a given site at T1 and T2
        # A single dissimilarity is computed because there are only two data vectors.
        #
        # Arguments --
        # vec1, vec2 : data vectors (species abundance or presence-absence data)
        # method == c("ruzicka", "%difference")
        # BCD=TRUE  : Compute and save the B and C components of the %difference and Ruzicka D.
        #             For the %difference, they are B/(2A+B+C), C/(2A+B+C), D/(2A+B+C).
        #             For the Ruzicka D, they are B/(A+B+C), C/(A+B+C), D/(A+B+C).
        # BCD=FALSE : Do not compute the components. BCD=FALSE for D other than %diff and Ruzicka.
        # ref=TRUE  : Compute the reference values of D, B and C
        #    =FALSE : Under permutation, compute only the value of D. Use separate code (shorter).
        #
        
        A <- sum(pmin(vec1, vec2))          # A = sum of minima from comparison of the 2 vectors
        sum.Y <- sum(vec1, vec2)            # Sum of all values in the two vectors, (2*A+B+C)
        #
        if (ref) {    # Compute the reference values of statistics D, B and C
            tmp <- vec1 - vec2
            B <- sum(tmp[tmp > 0])                 # Sum of the species losses between T1 and T2
            C <- -sum(tmp[tmp < 0])                # Sum of the species gains between T1 and T2
            D <- B + C                             # Dissimilarity
            
            # Under permutation, compute D from sum.Y and A; no need of the B and C values.
        } else { 
            D <- sum.Y - 2 * A                      # (2*A+B+C)-2*A = (B+C)
        }
        # Compute the denominator (den) of the Ruzicka or %difference index
        if (method == "ruzicka") 
            den <- (sum.Y - A)  # den = (A+B+C)
        else 
            den <- sum.Y                # den = (2A+B+C)
        if (!BCD) 
            B <- C <- NA
        list(B.den = B/den, C.den = C/den, D = D/den, B = B, C = C)
    }
    
    dissim <- function(mat1, mat2, n, method, BCD, ref) {
        # BCD=TRUE : Method is {"ruzicka", "%difference"} and output table BCD was requested
        # ref=TRUE : The function is called to compute the reference values of the TBI dissimil.
        
        vecD = vector(mode = "numeric", length = n)     # to receive the values D=(B+C)/den
        if (ref & BCD) { 
            vecB <- vector(mode = "numeric", length = n) # to receive the values B/den
            vecC <- vector(mode = "numeric", length = n) # to receive the values C/den
            v.B  <- vector(mode = "numeric", length = n) # to receive the values B
            v.C  <- vector(mode = "numeric", length = n) # to receive the values C
        } else {vecB <- vecC <- v.B <- v.C <- NA}
        #
        # Compute the dissimilarity between T1 and T2 for each object (site)
        #
        # 1. If method is "euclidean", "chord"
        #    compute the Euclidean distance
        if (any(method == c("euclidean", "chord")))  
            for (i in 1:n) vecD[i] <- dist(rbind(mat1[i,], mat2[i,])) 
            #
            # 2. Compute the Ruzicka or %difference dissimilarity 
            # "ruzicka"          # Quantitative form of Jaccard
            # "%difference"      # Quantitative form of Sorensen
            if (any(method == c("ruzicka", "%difference"))) { 
                for (i in 1:n) {
                    tmp <- RuzickaD(mat1[i,], mat2[i,], method = method, BCD = BCD, ref = ref) 
                    if (ref & BCD) {
                        vecB[i] <- tmp$B.den
                        vecC[i] <- tmp$C.den
                        v.B[i]  <- tmp$B
                        v.C[i]  <- tmp$C
                    }
                    vecD[i] <- tmp$D
                }
            }
            list(vecB = vecB, vecC = vecC, vecD = vecD, v.B = v.B, v.C = v.C)
    }
    
    ### End internal functions --
    
    A <- system.time({
        
        # Set "seed." to a specified integer, e.g. 1234, in function call to repeat a calculation
        if (!is.null(seed.)) set.seed(seed.)   
        
        epsilon <- sqrt(.Machine$double.eps)
        method <- match.arg(method, c("%difference", "ruzicka", "chord", "hellinger", "log.chord", "jaccard", "sorensen", "ochiai", "euclidean")) 
        n = nrow(mat1)
        p = ncol(mat1)
        if ((nrow(mat2) != n) | (ncol(mat2) != p)) stop("The matrices are not of the same size.")
        
        if (method == "hellinger") {
            mat1 <- sqrt(mat1) # sqrt() transformation done only once, before permutations
            mat2 <- sqrt(mat2) # sqrt() transformation done only once, before permutations
            method <- "chord" }
        if (method == "log.chord") { 
            mat1 <- log1p(mat1)    # log1p() transformation only once, before permutations
            mat2 <- log1p(mat2)    # log1p() transformation only once, before permutations
            method <- "chord" }
        if (method == "jaccard") {
            pa.tr <- TRUE
            method <- "ruzicka"
        }
        if (method == "sorensen") {
            pa.tr <- TRUE
            method <- "%difference"
        }
        if (method == "ochiai") {
            pa.tr <- TRUE
            method <- "chord"
        }
        #
        if (pa.tr) {
            mat1 <- ifelse(mat1 > 0, 1, 0)
            mat2 <- ifelse(mat2 > 0, 1, 0)
        }
        if (method == "chord") {
            tr <- TRUE
        } else { 
            tr <- FALSE
        }
        test.B.C <- NA 
        if ((any(method == c("ruzicka", "%difference"))) & BCD) { 
            BCD.mat <- matrix(0,n,3)
            if (method == "%difference")
                colnames(BCD.mat) <- c("B/(2A+B+C)","C/(2A+B+C)","D=(B+C)/(2A+B+C)")
            if (method == "ruzicka")    
                colnames(BCD.mat) <- c("B/(A+B+C)","C/(A+B+C)","D=(B+C)/(A+B+C)")
            rownames(BCD.mat) <- paste("Site", 1:n, sep = ".")
            Change <- vector(mode = "character", length = n)
        } else {
            BCD <- FALSE 
            BCD.mat <- NA 
            BCD.summ <- NA 
        }
        ###
        # 1. Compute the reference D for each object from corresponding vectors in the 2 matrices.
        if (tr) 
            tmp <- dissim(decostand(mat1, "norm"), decostand(mat2, "norm"), n, method, BCD, ref = TRUE)
        else 
            tmp <- dissim(mat1, mat2, n, method, BCD, ref = TRUE) 
        
        vecD.ref <- tmp$vecD
        BC <- NA
        if (BCD) { 
            BCD.mat[,1] <- tmp$vecB
            BCD.mat[,2] <- tmp$vecC
            BCD.mat[,3] <- tmp$vecD 
            for (i in 1:n) {
                if (tmp$vecB[i] > tmp$vecC[i]) Change[i] <- "-  " 
                else if (tmp$vecB[i] < tmp$vecC[i]) Change[i] <- "+  "
                else Change[i] <- "0  " 
            }
            BCD.summ <- matrix(NA, 1, 6)
            colnames(BCD.summ) <- c("mean(B/den)","mean(C/den)","mean(D)","B/(B+C)","C/(B+C)", 
                "Change")
            BCD.means <- apply(BCD.mat, 2, mean, na.rm = TRUE)  # Exclude the sites with value = NA
            BCD.summ[1,1:3] <- BCD.means
            BCD.summ[1,4:5] <- BCD.means[1:2] / BCD.means[3]
            BCD.summ <- as.data.frame(BCD.summ)
            if (BCD.summ[1,1] > BCD.summ[1,2]) BCD.summ[1,6] <- "-  " 
            else if (BCD.summ[1,1] < BCD.summ[1,2]) BCD.summ[1,6] <- "+  "
            else BCD.summ[1,6] <- "0  "
            rownames(BCD.summ) <- ""
            
            BCD.mat <- as.data.frame(BCD.mat)
            BCD.mat <- cbind(BCD.mat,Change)
            
            if ((n > 4) & test.BC) {   # Tests of significance of difference between B/den and C/den
                test.B.C <- matrix(NA,1,4)
                rownames(test.B.C) <- "Paired t.test"
                
                # Paired t-test between the vectors of B and C values
                if (test.t.perm) {
                    if (nperm < 999) nperm1 <- 999 else nperm1 <- nperm
                    t.res <- tpaired.randtest(tmp$vecC,tmp$vecB, nperm = nperm1, alternative = "two.sided", silent = TRUE)
                    test.B.C[1,] <- c(t.res$estim, t.res$t.ref, t.res$p.param, t.res$p.perm)
                    p.value <- t.res$p.param
                } else {
                    t.res <- t.test(tmp$vecC, tmp$vecB, paired = TRUE, alternative = "two.sided")
                    test.B.C[1,] <- c(t.res$estimate, t.res$statistic, t.res$p.value, NA)
                    p.value <- t.res$p.value
                }
                signif. <- ifelse(p.value > 0.05, " ", "*")
                test.B.C <- as.data.frame(test.B.C)
                test.B.C <- cbind(test.B.C, signif.)
                colnames(test.B.C) <- c("  mean(C-B)", "Stat", "p.param", "p.perm","  p<=0.05")
            } else {
                test.B.C <- NA
            }
            # Matrix containing the observed values of B and C, in case they are needed later
            if (save.BC) {
                BC <- cbind(tmp$v.B, tmp$v.C)
                colnames(BC) <- c("B", "C")
                rownames(BC) <- paste("Site", 1:n, sep = ".")
            } 
        }
        
        ###
        # 2. Permutation method: permute data separately in each column. 
        # Permute *the raw data* by columns. Permute the two matrices in the same way, 
        # saving the seed before the two sets of permutations through sample(). 
        # There is a separate permutation test for each distance in vector D.
        if (nperm > 0) {
            my.vec <- sample(1:(10 * nperm), size = nperm)
            # cat("seed =",my.vec,'\n')
            BCD <- FALSE
            nGE.D = rep(1,n)
            for (iperm in 1:nperm) {			 
                set.seed(my.vec[iperm])
                mat1.perm <- apply(mat1, 2, sample, replace = replace)
                set.seed(my.vec[iperm])
                mat2.perm <- apply(mat2, 2, sample, replace = replace)
                # 3. Recompute transformations of the matrices 
                #    and the dissimilarity value of the paired vectors of each site i.
                if (tr) {
                    tmp <- dissim(decostand(mat1.perm, "norm"), 
                        decostand(mat2.perm, "norm"), n, method, BCD, ref = FALSE)
                } else {
                    tmp <- dissim(mat1.perm, mat2.perm, n, method, BCD, ref = FALSE)
                }
                vecD.perm <- tmp$vecD
                ge <- which(vecD.perm + epsilon >= vecD.ref)
                if (length(ge) > 0) nGE.D[ge] <- nGE.D[ge] + 1
            }
            # 4. Compute the p-value associated with each distance (i.e. site).
            p.dist <- nGE.D / (nperm + 1)
        } else
            p.dist <- NA   # if nperm=0
        
        p.adj <- p.adjust(p.dist,"holm")
    })
    A[3] <- sprintf("%2f",A[3])
    if (clock) cat("Computation time =", A[3]," sec",'\n')
    #
    out <- list(TBI = vecD.ref, p.TBI = p.dist, p.adj = p.adj, BCD.mat = BCD.mat, BCD.summary = BCD.summ, t.test_B.C = test.B.C, BC = BC)
    class(out) <- "TBI"
    return(out)
}
