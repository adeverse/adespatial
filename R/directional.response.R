#' Directional indices of community change
#'
#' Compute directional indices of community change along coenoclines or time.
#'
#' @param mat A community composition data matrix with sites in rows and species 
#' in columns. The direction of the physical process is indicated by the order of the 
#' sampling units in \code{mat}. The class of \code{mat} can be either \code{data.frame} 
#' or \code{matrix}.
#' 
#' @param method One of the 11 calculation methods available in the function:
#' \code{"overlap"}, \code{"gain"}, \code{"loss"}, \code{"gaining.turnover"}, 
#' \code{"neutral.turnover"}, \code{"losing.turnover"}, \code{"gaining.nestedness"}, 
#' \code{"neutral.nestedness"}, \code{ "losing.nestedness"}, 
#' \code{"gaining.strict.nestedness"}, \code{"losing.strict.nestedness"}. 
#' The default value is \code{method="overlap"}.
#' 
#' @param relativize Compute relativized indices: \code{relativize="J"} for the
#' Jaccard denominator (a+b+c) or \code{relativize="S"} for the Sorensen denominator 
#' (2*a+b+c). If \code{relativize=NULL}, the index is not divided by a denominator.
#' 
#' @details
#' The output matrix is non-symmetric, meaning that its upper triangle is not the mirror 
#' image of the lower triangle. Given the direction of the physical process through space 
#' or time indicated by the order of the sampling units, the output matrix \code{mat.out} 
#' reflects that direction in its non-symmetric presentation, \emph{From} (rows of the 
#' matrix) \emph{To} (columns of the matrix). Users of the function can extract one or the 
#' other of these triangular portions and analyse them separately. See Appendix xx for 
#' examples.
#'
#' @return A list containing the following results: \itemize{ \item
#' \code{mat.out}: A square matrix with the chosen index computed among the sites.
#'	 Depending on the method chosen, this matrix may be symmetric or non-symmetric.
#' \code{total.t}: methods  #4 to 6, a matrix with total turnover (b+c); else NA.
#' \code{total.n}: For methods #7 to 9, a matrix with total nestedness a+abs(b-c) 
#'   if a>0; else NA.
#' \code{total.strict.n}: For methods #10 and 11, a matrix with total strict nestedness 
#'   a+abs(b-c) if a>0 and b!=c; else NA.
#' \code{den}: For calculation results with Jaccard or Sorensen denominators: 
#' 	 a square matrix of denominators. If \code{relativize=NULL}, \code{den=NA}. }
#' 
#' @author Dénes Schmera \email{schmera.denes@blki.hu} 
#' and Pierre Legendre \email{pierre.legendre@umontreal.ca}
#'
#' @references
#' Schmera, D., P. Legendre, T. Eros, M. Toth, E. K. Magyari, B. Baur and J. Podani. 
#' 2022. New measures for quantifying directional changes in presence-absence community 
#' data. Ecological Indicators 136: 108618. https://doi.org/10.1016/j.ecolind.2022.108618
#'
#' Verneaux, J. (1973) \emph{Cours d'eau de Franche-Comté (Massif du Jura). 
#' Recherches écologiques sur le réseau hydrographique du Doubs. Essai de biotypologie}. 
#' Thèse d'État, Besançon. 1–257.
#' 
#' @examples
#'
#' # Artificial Example
#' art <- c(1,1,1,0,0,0,
#'          0,0,0,1,1,0,
#'          0,0,0,0,0,1)
#' art.data <- matrix(art, nrow=3, ncol=6, byrow=TRUE)
#'
#' art.out <- directional.response(art.data, method="overlap",relativize=NULL)
#' 
#' # Real data example: the Doubs River fish data (Verneaux 1973), available in ade4.
#' # 30 sites, 27 species. No fish had been caught at site 8; remove that site
#' if(require("ade4", quietly = TRUE)) {
#' 
#' data(doubs)
#' dim(doubs$fish)   
#' fish <- doubs$fish[-8,] 
#' dim(fish)
#' doubs.out <- directional.response(fish, method="gain", relativize="S")
#' }
#' 
#' @export directional.response
#'


directional.response <-
	function(mat, 
	method="overlap", 
	relativize=NULL)
{
METHODS <- c("overlap","gain","loss", "gaining.turnover","neutral.turnover", "losing.turnover", "gaining.nestedness","neutral.nestedness", "losing.nestedness", "gaining.strict.nestedness", "losing.strict.nestedness")
#
method <- pmatch(method, METHODS)   # Rank order of the selected method in list METHODS
inm <- METHODS[method]              # Name of the selected method
if (is.na(method)) 
        stop("Invalid method")

mat <- as.matrix(mat)
mat.b <- ifelse(mat > 0, 1, 0)  # Make the data matrix binary, i.e. 0/1 data

a <- mat.b %*% t(mat.b)         # overlap
c <- mat.b %*% (1 - t(mat.b))   # gain   # b in beta.div.comp.R
b <- (1 - mat.b) %*% t(mat.b)   # loss   # c in beta.div.comp.R
min.bc <- pmin(b,c)
repl <- 2*min.bc
rich <- abs(b-c)
change <- (b-c)   # positive values for species losses, negative for gains

zeromatrix<-matrix(rep(0,nrow(mat.b)*ncol(mat.b)),nrow=nrow(mat.b),ncol=ncol(mat.b))
ga<-ifelse(change<0,rich, zeromatrix)
lo<-ifelse(change>0,rich, zeromatrix)
gt<-repl+ga
nt<-repl
lt<-repl+lo
gn<-ifelse(a==0,zeromatrix,ifelse(change>=0,a,a+rich))
nn<-ifelse(a==0,zeromatrix,a)
ln<-ifelse(a==0,zeromatrix,ifelse(change<=0,a,a+rich))
gsn<-ifelse(a==0|change==0,zeromatrix,ifelse(change>0,a,a+rich))
lsn<-ifelse(a==0|change==0,zeromatrix,ifelse(change<0,a,a+rich))

if (method == 1) out1<-a     # overlap
if (method == 2) out1<-c     # gain
if (method == 3) out1<-b     # loss
if (method == 4) out1<-gt    # gt
if (method == 5) out1<-nt    # nt
if (method == 6) out1<-lt    # lt
if (method == 7) out1<-gn    # gn
if (method == 8) out1<-nn    # nn
if (method == 9) out1<-ln    # ln
if (method == 10) out1<-gsn  # gsn
if (method == 11) out1<-lsn  # lsn
# Compute total.t, noted tt, for comparison with outpput of methods 4 to 6
if(method>3 & method<7) tt <- b+c else tt <- NA
# Compute total.n, noted tn, for for comparison with output of methods 7 to 9
if(method>6 & method<10) tn <- ifelse(a>0,a+abs(b-c),0) else tn <- NA
# Compute total.strict.n, noted tn2, for for comparison with output of methods 10 to 11
if(method>9) tn2 <- ifelse(a>0,ifelse(b==c,0,a+abs(b-c)),0) else tn2 <- NA

if (is.null(relativize)) { 
	out2 <- out1 
	den <- NA
	} else {
		if(relativize == "J") {
			den <- a+b+c
			out2 <- out1/den
		}
		if(relativize == "S") {
			den <- 2*a+b+c
			#
			if(method == 1) {          # overlap
				out2 <- 2*a/den
			} else if(method == 2) {   # gain
				out2 <- c/den
			} else if(method == 3) {   # loss
				out2 <- b/den
			} else if(method == 4) {   # gt
				out2 <- gt/den
			} else if(method == 5) {   # nt
				out2 <- nt/den	
			} else if(method == 6) {   # lt
				out2 <- lt/den
			} else if(method == 7) {   # gn
				out2 <- ifelse(a==0,zeromatrix,ifelse(change>=0,2*a/den,(2*a+rich)/den))
			} else if(method == 8) {   # nn
				out2 <- ifelse(a==0,zeromatrix,2*a/den)
			} else if(method == 9) {   # ln
				out2 <- ifelse(a==0,zeromatrix,ifelse(change<=0,2*a/den,(2*a+rich)/den))
			} else if(method == 10) {  # gsn
				out2 <- ifelse(a==0|change==0,zeromatrix,
				ifelse(change>0,2*a/den,(2*a+rich)/den))
			} else if(method == 11) {  # lsn
				out2 <- ifelse(a==0|change==0,zeromatrix,
				ifelse(change<0,2*a/den,(2*a+rich)/den))
			}
		} 
}
#			
cat("Method: ",inm,"\n")
#
if(!is.null(relativize)) {
	if(relativize == "J") {
		cat("Relativize with denominator: Jaccard\n") 
	} else if(relativize == "S") {
		cat("Relativize with denominator: Sorensen\n") 
	}
}
#
list(mat.out=out2, total.t=tt, total.n=tn, total.strict.n=tn2, den=den)
}
