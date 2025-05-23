% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dist.ldc.R
\name{dist.ldc}
\alias{dist.ldc}
\title{Dissimilarity matrices for community composition data}
\usage{
dist.ldc(Y, method = "hellinger", binary = FALSE, samp = TRUE, silent = FALSE)
}
\arguments{
\item{Y}{Community composition data. The object class can be either
\code{data.frame} or \code{matrix}.}

\item{method}{One of the 21 dissimilarity coefficients available in the
function: \code{"hellinger"}, \code{"chord"}, \code{"log.chord"},
\code{"chisquare"}, \code{"profiles"}, \code{"percentdiff"},
\code{"ruzicka"}, \code{"divergence"}, \code{"canberra"},
\code{"whittaker"}, \code{"wishart"}, \code{"kulczynski"},
\code{"jaccard"}, \code{"sorensen"}, \code{"ochiai"}, \code{"ab.jaccard"},
\code{"ab.sorensen"}, \code{"ab.ochiai"}, \code{"ab.simpson"},
\code{"euclidean"}, \code{"manhattan"}, \code{"modmeanchardiff"}. See
Details. Names can be abbreviated to a non-ambiguous set of first letters.
Default: \code{method="hellinger"}.}

\item{binary}{If \code{binary=TRUE}, the data are transformed to
presence-absence form before computation of the dissimilarities. Default
value: \code{binary=FALSE}, except for the Jaccard, Sørensen and Ochiai
indices where \code{binary=TRUE}.}

\item{samp}{If \code{samp=TRUE}, the abundance-based distances (ab.jaccard,
ab.sorensen, ab.ochiai, ab.simpson) are computed for sample data. If
\code{samp=FALSE}, binary indices are computed for true population data.}

\item{silent}{If \code{silent=FALSE}, informative messages sent to users will
be printed to the R console. Use \code{silent=TRUE} is called on a
numerical simulation loop, for example.}
}
\value{
A dissimilarity matrix, with class \code{dist}.
}
\description{
Compute dissimilarity indices for ecological data matrices. The dissimilarity
indices computed by this function are those described in Legendre and De
Cáceres (2013). In the name of the function, 'ldc' stands for the author's
names. Twelve of these 21 indices are not readily available in other R
package functions; four of them can, however, be computed in two computation
steps in \code{vegan}.
}
\details{
The dissimilarities computed by this function are the following.
  Indices i and k designate two rows (sites) of matrix Y, j designates a
  column (species). D[ik] is the dissimilarity between rows i and k. p is the
  number of columns (species) in Y; pp is the number of species present in
  one or the other site, or in both. y[i+] is the sum of values in row i;
  same for y[k+]. y[+j] is the sum of values in column j. y[++] is the total
  sum of values in Y. The indices are computed by functions written in C for
  greater computation speed with large data matrices. \itemize{ 
  \item Group 1 - D computed by transformation of Y followed by Euclidean distance
  \itemize{ \item Hellinger D, D[ik] =
  sqrt(sum((sqrt(y[ij]/y[i+])-sqrt(y[kj]/y[k+]))^2)) (Rao 1995)
  \item chord D, D[ik] =
  sqrt(sum((y[ij]/sqrt(sum(y[ij]^2))-y[kj]/sqrt(sum(y[kj]^2)))^2)) 
  (Orlóci 1967)
  \item log-chord D, D[ik] = chord D[ik] computed on 
  log(y[ij]+1)-transformed data (Legendre & Borcard 2018) 
  \item chi-square D, D[ik] = sqrt(y[++]  
  sum((1/j[+j])(y[ij]/y[i+]-y[kj]/y[k+])^2)) (Lebart & Fénelon 1971)
  \item species profiles D, D[ik] = sqrt(sum((y[ij]/y[i+]-y[kj]/y[k+])^2)) 
  (Legendre & Legendre (2012) }

  \item Group 2 - Other D functions appropriate for beta diversity studies
  where A = sum(min(y[ij],y[kj])), B = y[i+]-A, C = y[k+]-A 
  \itemize{ \item percentage difference D (aka Bray-Curtis), 
  D[ik] = (sum(abs(y[ij]-y[k,j])))/(y[i+]+y[k+]) 
  or else, D[ik] = (B+C)/(2A+B+C) (Odum 1950)
  \item Ružička D, D[ik] = 1-(sum(min(y[ij],y[kj])/sum(max(y[ij],y[kj]))
  or else, D[ik] = (B+C)/(A+B+C) (Ružička 1958)
  \item coefficient of divergence D, D[ik] =
  sqrt((1/pp)sum(((y[ij]-y(kj])/(y[ij]+y(kj]))^2)) (Clark 1952)
  \item Canberra metric D,
  D[ik] = (1/pp)sum(abs(y[ij]-y(kj])/(y[ij]+y(kj])) (Lance & Williams 1967)
  \item Whittaker D, 
  D[ik] = 0.5*sum(abs(y[ij]/y[i+]-y(kj]/y[k+])) (Whittaker 1952)
  \item Wishart D, D[ik] = 1-sum(y[ij]y[kj])
  /(sum(y[ij]^2)+sum(y[kj]^2)-sum(y[ij]y[kj])) (Wishart 1969)
  \item Kulczynski D, D[ik] =
  1-0.5((sum(min(y[ij],y[kj])/y[i+]+sum(min(y[ij],y[kj])/y[k+])) }

  \item Group 3 - Classical indices for binary data; they are appropriate for
  beta diversity studies. Value a is the number of species found in both i
  and k, b is the number of species in site i not found in k, and c is the
  number of species found in site k but not in i. The D matrices are
  square-root transformed, as in dist.binary of ade4; the user-oriented
  reason for this transformation is explained below. \itemize{
  \item Jaccard coefficient of community (CC), D[ik] = sqrt((b+c)/(a+b+c))
  (Jaccard 1901, CC described in a sentence; 1908, CC shown as an equation)
  \item Sørensen D, D[ik] = sqrt((b+c)/(2a+b+c)) (Sørensen 1948)
  \item Ochiai D, D[ik] = sqrt(1 - a/sqrt((a+b)(a+c))) (Ochiai 1957) }

  \item Group 4 - Abundance-based indices of Chao et al. (2006) for
  quantitative abundance data. These functions correct the index for species
  that have not been observed due to sampling errors. For the meaning of the
  U and V notations, see Chao et al. (2006, section 3). When
  \code{samp=TRUE}, the abundance-based distances (ab.jaccard, ab.sorensen,
  ab.ochiai, ab.simpson) are computed for sample data. If \code{samp=FALSE},
  indices are computed for true population data. - Do not use indices of
  group 4 with \code{samp=TRUE} on presence-absence data; the indices are not
  meant to accommodate this type of data. If \code{samp=FALSE} is used with
  presence-absence data, the indices are the regular Jaccard, Sørensen,
  Ochiai and Simpson indices. On output, however, the D matrices are not
  square-rooted, contrary to the Jaccard, Sørensen and Ochiai indices in
  section 3 which are square-rooted. \itemize{ 
  \item abundance-based Jaccard D, D[ik] = 1-(UV/(U+V-UV)) 
  \item abundance-based Sørensen D, D[ik] = 1-(2UV/(U+V)) 
  \item abundance-based Ochiai D, D[ik] = 1-sqrt(UV) 
  \item abundance-based Simpson D, D[ik] = 1-(UV/(UV+min((U-UV),(V-UV)))) }

  \item Group 5 - General-purpose dissimilarities that do not have an upper
  bound (maximum D value). They are inappropriate for beta diversity studies.
  \itemize{ 
  \item Euclidean D, D[ik] = sqrt(sum(y[ij]-y[kj])^2) 
  \item Manhattan D, D[ik] = sum(abs(y[ij] - y[ik])) 
  \item modified mean character difference, 
  D[ik] = (1/pp) sum(abs(y[ij] - y[ik])) (Legendre & Legendre (2012) } } 

  The properties of
  all dissimilarities available in this function (except Ružička D) were
  described and compared in Legendre & De Cáceres (2013), who showed that
  most of these dissimilarities are appropriate for beta diversity studies.
  Inappropriate are the Euclidean, Manhattan, modified mean character
  difference, species profile and chi-square distances. Most of these
  dissimilarities have a maximum value of either 1 or sqrt(2). Three
  dissimilarities (Euclidean, Manhattan, Modified mean character difference)
  do not have an upper bound and are thus inappropriate for beta diversity
  studies. The chi-square distance has an upper bound of sqrt(2*(sum(Y))).
  \cr\cr The Euclidean, Hellinger, chord, chi-square and
  species profiles dissimilarities have the property of being Euclidean,
  meaning that they never produce negative eigenvalues in principal
  coordinate analysis. The Canberra, Whittaker, percentage difference,
  Wishart and Manhattan coefficients are Euclidean when they are square-root
  transformed (Legendre & De Cáceres 2013, Table 2). The distance forms (1-S)
  of the Jaccard, Sørensen and Ochiai similarity (S) coefficients are
  Euclidean after taking the square root of (1-S) (Legendre & Legendre 2012,
  Table 7.2). The D matrices resulting from these three coefficients are
  outputted in the form sqrt(1-S), as in function \code{dist.binary} of ade4,
  because that form is Euclidean and will thus produce no negative
  eigenvalues in principal coordinate analysis. \cr\cr The Hellinger, chord,
  chi-square and species profile dissimilarities are computed using the
  two-step procedure developed by Legendre & Gallagher (2001). The data are
  first transformed using either the row marginals, or the row and column
  marginals in the case of the chi-square distance. The dissimilarities are
  then computed from the transformed data using the Euclidean distance
  formula. As a consequence, these four dissimilarities are necessarily
  Euclidean. D matrices for other binary coefficients can be computed in two
  ways: either by using function \code{dist.binary} of ade4, or by choosing
  option \code{binary=TRUE}, which transforms the abundance data to binary
  form, and using one of the quantitative indices of the present function.
  Table 1 of Legendre & De Cáceres (2013) shows the incidence-based
  (presence-absence-based) indices computed by the various indices using
  binary data. \cr\cr The Euclidean distance computed on untransformed
  presence-absence or abundance data produces non-informative and incorrect
  ordinations, as shown in Legendre & Legendre (2012, p. 300) and in Legendre
  & De Cáceres (2013). However, the Euclidean distance computed on
  log-transformed abundance data produces meaningful ordinations in principal
  coordinate analysis (PCoA). Nonetheless, it is easier to compute a PCA of
  log-transformed abundance data instead of a PCoA; the resulting ordination
  with scaling 1 will be meaningful. Messages are printed to the R console
  indicating the Euclidean status of the computed dissimilarity matrices.
  Note that for the chi-square distance, the columns that sum to zero are
  eliminated before calculation of the distances, thus preventing divisions
  by zero in the calculation of the chi-square transformation.
}
\examples{

if(require("vegan", quietly = TRUE)) {
data(mite)
mat1  = as.matrix(mite[1:10, 1:15])   # No column has a sum of 0
mat2 = as.matrix(mite[61:70, 1:15])   # 7 of the 15 columns have a sum of 0

#Example 1: compute Hellinger distance for mat1
D.out = dist.ldc(mat1,"hellinger")

#Example 2: compute chi-square distance for mat2
D.out = dist.ldc(mat2,"chisquare")

#Example 3: compute percentage difference dissimilarity for mat2
D.out = dist.ldc(mat2,"percentdiff")

}


}
\references{
Chao, A., R.L. Chazdon, R.K. Colwell and T.J. Shen. 2006.
  Abundance-based similarity indices and their estimation when there are
  unseen species in samples. Biometrics 62: 361-371.

  Clark, P. J. 1952. An extension of the coefficient of divergence for 
  use with multiple characters. Copeia 1952: 61-64.

  Jaccard, P. 1901. Étude comparative de la distribution florale dans  
  dune portiones Alpes et du Jura. Bulletin de la Société vaudoise des sciences 
  naturelles 37: 547-579.

  Jaccard, P. 1908. Nouvelles recherches sur la distribution florale. Bulletin 
  de la Société vaudoise des sciences naturelles 44: 223-270.

  Kulczynski, S. 1928. Die Pflanzenassoziationen der Pieninen. Bull. Int. 
  Acad. Pol. Sci. Lett. Cl. Sci. Math. Nat. Ser. B, Suppl. II (1927): 57–203. 

  Lance, G. N. and W. T. Williams. 1967. Mixed-data classificatory programs. I. 
  Agglomerative systems. Australian Computer Journal 1: 15-20. 

  Lebart, L. and J. P. Fénelon. 1971. Statistique et informatique appliquées. Dunod, Paris.

  Legendre, P. and D. Borcard. 2018. Box-Cox-chord transformations for
  community composition data prior to beta diversity analysis. Ecography 41:
  1820-1824.

  Legendre, P. and M. De Cáceres. 2013. Beta diversity as the variance of
  community data: dissimilarity coefficients and partitioning. Ecology
  Letters 16: 951-963.

  Legendre, P. and E. D. Gallagher, E.D. 2001. Ecologically meaningful
  transformations for ordination of species data. Oecologia 129: 271-280.

  Legendre, P. and Legendre, L. 2012. Numerical Ecology, 3rd English edition.
  Elsevier Science BV, Amsterdam.

  Ochiai, A. 1957. Zoogeographic studies on the soleoid fishes found in Japan 
  and its neighbouring regions–II. Bulletin of the Japanese Society of Scientific 
  Fisheries 22: 526-530.

  Odum, E. P. 1950. Bird populations of the Highlands (North Carolina) plateau 
  in relation to plant succession and avian invasion. Ecology 31: 587-605. 

  Orlóci, L. 1967. An agglomerative method for classification of plant communities. 
  Journal of Ecology 55: 193-206.

  Rao, C. R. 1995. A review of canonical coordinates and an alternative to 
  correspondence analysis using Hellinger distance. Qüestiió (Quaderns 
  d’Estadística i Investigació Operativa) 19: 23-63.

  Ružička, Milan. 1958. Anwendung mathematisch-statistiker Methoden in Geobotanik 
  (synthetische Bearbeitung von Aufnahmen). Biologia (Bratislava) 13: 647–661. 

  Sørensen, T. 1948. A method of establishing groups of equal amplitude in plant 
  sociology based on similarity of species content and its application to analysis 
  of the vegetation on Danish commons. Biologiske Skrifter 5: 1-34.

  Whittaker, R. H. 1952. A study of summer foliage insect communities in the 
  Great Smoky Mountains. Ecological Monographs 22: 1–44.

  Wishart, D. (1969). CLUSTAN 1a User Manual. Computing Laboratory, University 
  of St. Andrews, St. Andrews, Fife, Scotland.
}
\author{
Pierre Legendre \email{pierre.legendre@umontreal.ca} and Naima Madi
}
