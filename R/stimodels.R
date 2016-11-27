#' Space-time interaction in ANOVA without replication
#' 
#' Function \code{stimodels} performs two-way ANOVA to test space-time
#' interaction without replicates using one among a set of possible models.
#' Function \code{quicksti} allows performing space-time ANOVA in a simplified
#' way. In many models, degrees of freedom are saved by coding space and/or
#' time parsimoniously using distance-based Moran Eigenvector Maps (dbMEM).
#' 
#' @details In \code{stimodels} tests for space-time interaction and space or time main
#' effects are conducted using one of the different models. With Models 2, 6a
#' and 6b the interaction test is not available.
#' 
#' Model 2 - Space and Time are coded using Helmert contrasts for the main
#' effects. No interaction is tested. Model 3a - Space is coded using dbMEM
#' variables whereas Time is coded using Helmert contrasts. Model 3b - Space is
#' coded using Helmert contrasts whereas Time is coded using dbMEM variables.
#' Model 4 - Both Space and Time are coded using dbMEM variables for all tests.
#' Model 5 - Space and Time are coded using Helmert contrasts for the main
#' factor effects, but they are coded using dbMEM variables for the interaction
#' term. Model 6a - Nested model. Testing for the existence of spatial
#' structure (common or separate) using dbMEM variables to code for Space.
#' Model 6b - Nested model. Testing for the existence of temporal structure
#' (common or separate) using dbMEM variables to code for Time. Model 7 -
#' Space and Time are coded using dbMEM variables for the main factor effects,
#' but they are coded using Helmert contrasts for the interaction term (not
#' recommended).
#' 
#' When using \code{quicksti}, space-time interaction is first tested using
#' Model 5. Depending on the outcome of this test, the main factors are tested
#' using different strategies. If the interaction is not significant then the
#' test of main factors is also done following Model 5. If the interaction is
#' significant, then a nested model (6a) is used to know whether separate
#' spatial structures exist and another (6b) to know whether separate temporal
#' structures exist. In \code{quicksti} function space and time are always
#' considered fixed factors (F ratios are constructed using residual MS in the
#' denominator).
#' 
#' For the interaction the permutations are unrestricted, whereas for the main
#' factors the permutations are restricted within time blocks (for the test of
#' factor Space) or space blocks (for the test of factor Time). By default, the
#' function computes dbMEM for space and time coding, but other space and/or
#' time descriptors can be provided by the user, through \code{COD.S} and
#' \code{COD.T}.
#' 
#' 
#' @param Y Site-by-species response data table. Assumes row blocks
#' corresponding to times, i.e. within each block all sites are provided 
#' (in the same order).
#' @param S Number of spatial points or a matrix of spatial coordinates.
#' @param Ti Number of time campaigns or a matrix (a vector) of temporal
#' coordinates.
#' @param model Linear space-time model to be used (can be either "2", "3a",
#' "3b", "4", "5", "6a", "6b", or "7").
#' @param nperm Number of permutations in the significance tests.
#' @param alpha In \code{quicksti}, confidence level for the interaction test. 
#' Depending on the decision for the interaction test, the main factors are 
#' tested differently.
#' @param nS Number of space dbMEMs to use (by default, -1, all dbMEMs with
#' positive autocorrelation are used).
#' @param nT Number of time dbMEMs to use (by default, -1, all dbMEMs with
#' positive autocorrelation are used).
#' @param Sfixed Logical: is factor Space fixed, or not (if FALSE, it is
#' considered a random factor).
#' @param Tfixed Logical: is factor Time fixed, or not (if FALSE, it is
#' considered a random factor).
#' @param COD.S Spatial coding functions to be used instead of dbMEM. The
#' number of columns must be lower than \code{S} and the number of rows equal
#' to the number of rows in \code{Y}.
#' @param COD.T Temporal coding functions to be used instead of dbMEM. The
#' number of columns must be lower than \code{Ti} and the number of rows equal
#' to the number of rows in \code{Y}.
#' @param print.res If TRUE displays the results and additional information
#' onscreen (recommended).
#' @return
#' 
#' \item{testS}{ An object with the result of the space effect test, including
#' the mean squares for the F numerator (\code{MS.num}), the mean squares for
#' the F denominator (\code{MS.den}), the proportion of explained variance
#' (\code{R2}), the adjusted proportion of explained variance (\code{R2.adj}),
#' the F statistics (\code{F}) and its p-value computed from a permutation test
#' (\code{Prob}). } \item{testT}{ An object with the result of the time effect
#' test, like \code{testS}.} \item{teststi}{ An object with the result of the
#' space-time interaction test, like \code{testS}.}
#' @author Pierre Legendre \email{pierre.legendre@@umontreal.ca}, Miquel De Caceres and Daniel Borcard
#' @seealso \code{\link{trichoptera}}
#' @references Borcard, D. and P. Legendre. 2002. All-scale spatial analysis of
#' ecological data by means of principal coordinates of neighbour matrices.
#' Ecological Modelling 153: 51-68.
#' 
#' Dray, S., P. Legendre and P. R. Peres-Neto. 2006. Spatial modelling: a
#' comprehensive framework for principal coordinate analysis of neighbour
#' matrices (PCNM). Ecological Modelling 196: 483-493.
#' 
#' Legendre, P., M. De Caceres and D. Borcard. 2010. Community surveys through
#' space and time to assess environmental changes: testing space-time
#' interaction in the absence of replication. Ecology 91: 262-272.
#' @keywords multivariate spatial
#' @examples
#' 
#' data(trichoptera)
#' 
#' # log-transform species data (excluding site and time colums)
#' trich.log <- log1p(trichoptera[,-c(1,2)]) 
#' 
#' 
#' # Run space-time interaction test using model "5"
#' stimodels(trich.log, S=22, Ti=10, nperm=99, model="5")
#' 
#' # Run space-time analysis with tests for main effects after testing 
#' # interaction (which is significant)
#' quicksti(trich.log, S=22, Ti=10, nperm=99)
#' 
#' # Run space-time analysis for time blocks number 6 and 7. 
#' # Interaction is then not significant and tests of main effects are done 
#' # following model 5
#' quicksti(trich.log[111:154,], S=22, Ti=2, nperm=999)
#' 
#' @importFrom MASS ginv
#' @importFrom stats model.matrix
#' @aliases stimodels quicksti
#' @export stimodels quicksti
#' @rdname stimodels  
#'   

'stimodels' <- function(Y, S, Ti,  model="5", nperm=999, nS=-1, nT=-1, Sfixed=TRUE, Tfixed=TRUE, COD.S=NULL, COD.T=NULL, print.res=TRUE)
{
	if(!is.logical(print.res)) {
		stop("Wrong operator; 'print.res' should be either 'FALSE' or 'TRUE'", call.=FALSE)
	}
	if(model!="2" && model!="3a" && model!="3b" && model!="4" && model!="4" && model!="5" && model!="6a" && model!="6b" && model!="7") {
		stop(paste("Unrecognized model ",model,"; 'model' should be '2', '3a', '3b', '4', '5', '6a','6b' or '7'.", sep=""), call.=FALSE)
	}
	
	aa <- system.time( {
	    # Sets the number and location of spatial and temporal points
		S <- as.matrix(S)
		if(dim(S)[1]==1 && dim(S)[2]==1) {
			s <- S[1,1]
			sitesX <- c(1:s)
		} else {
			s <- dim(S)[1]
			sitesX <- S
		}		
		Ti <- as.matrix(Ti)
		if(dim(Ti)[1]==1 && dim(Ti)[2]==1) {
			tt <- Ti[1,1]
			timesX <- c(1:tt) 
		} else {
			tt <- dim(Ti)[1]
			timesX <- Ti
		}
		
		# Total number of rows in matrix
		n <- s*tt			
		p <- dim(Y)[2]

		# Check response data file containing species data
		Y <- as.matrix(Y)
		p <- dim(Y)[2]
		
		if(dim(Y)[1] != n) stop("The number of rows in species file is not (S x Ti)", call.=FALSE) 

		# Center response data
		Y <- scale(Y, center=TRUE, scale=FALSE)

		
		if(print.res) {
			cat("=======================================================\n")
			cat("        Space-time ANOVA without replicates\n")
			cat("                                                  \n")
			cat("  Pierre Legendre, Miquel De Caceres, Daniel Borcard\n")
			cat("=======================================================\n\n")
			cat(" Number of space points (s) =", s,'\n')
			cat(" Number of time points (tt) =", tt,'\n')
			cat(" Number of observations (n = s*tt) =", n,'\n')
			cat(" Number of response variables (p) =", p,'\n','\n')
		}

		# Generates space dbMEM variables (if necessary)
		if(model=="3a"|| model=="6a" || model=="7"|| model=="4"|| model=="5") {
			if(is.null(COD.S)) {	# Generate spatial dbMEMs if not given by user
				if(print.res) cat(" Computing dbMEMs to code for space\n")
				dbMEM.S.tmp <- dbmem(sitesX, MEM.autocor="positive")
				SS <- as.matrix(dbMEM.S.tmp)
				dbMEM.S.thresh <- give.thresh(dist(sitesX))
				if(print.res) cat(" Truncation level for space dbMEMs =", dbMEM.S.thresh, "\n")
				dbMEM.S <- SS
				for(j in 2:tt) dbMEM.S <- rbind(dbMEM.S,SS)
				if(nS==-1) nS <- ncol(SS)
				else {
				   if(nS > ncol(SS)) {
				      cat("\n Number of requested spatial variables nS =", nS, "larger than available\n")
				      cat(" The", ncol(SS), "available spatial dbMEM will be used\n\n")
				      nS <- ncol(SS)
				      }
				    else dbMEM.S <- dbMEM.S[,1:nS]
				      }
			} else {
				dbMEM.S <- apply(as.matrix(COD.S),2,scale,center=TRUE,scale=TRUE)
				nS <- dim(dbMEM.S)[2]
				if(nS>=s) stop("The number of spatial coding functions must be lower than s", call.=FALSE) 
			  if(nrow(as.matrix(COD.S))!=dim(Y)[1]) stop("The number of rows in COD.S must be equal to the number of observations", call.=FALSE) 
			}
		}
		
		# Generate dbMEM variables for time (if necessary)
		if(model=="3b" || model=="6b"|| model=="7"||model=="4" || model=="5") {
			if(is.null(COD.T)) {      # Generate temporal dbMEM variables if not given by user
				if(print.res) cat(" Computing dbMEMs to code for time\n")
				dbMEM.T.tmp <- dbmem(timesX, MEM.autocor="positive")
				TT <- as.matrix(dbMEM.T.tmp)
				dbMEM.T.thresh <- give.thresh(dist(timesX))
				if(print.res) cat(" Truncation level for time dbMEMs =", dbMEM.T.thresh, "\n\n")
				T.temp <- TT[1,]
				for(i in 2:s) T.temp <- rbind(T.temp,TT[1,])
				dbMEM.T <- as.matrix(T.temp)
				for(j in 2:tt) {
					T.temp <- TT[j,]
					for(i in 2:s) T.temp <- rbind(T.temp,TT[j,])
					dbMEM.T <- as.matrix(rbind(dbMEM.T,T.temp))
				}
				if(nT==-1)  nT <- ncol(TT)
				else {
				   if(nT > ncol(TT)) {
				      cat("\n Number of requested temporal variables nT =", nT, "larger than available\n")
				      cat(" The", ncol(TT), "available temporal dbMEM will be used\n\n")
				      nT <- ncol(TT)
				      }
				    else dbMEM.T <- dbMEM.T[,1:nT]
				      }
				if(nT) {
				   if(nT > ncol(TT)) {
				      cat("Number of requested temporal variables nT=", nT, "larger than available\n")
				      cat("The", ncol(TT), "available temporal dbMEM will be used\n\n")
				      nT <- ncol(TT)
				      }
				    else dbMEM.T <- dbMEM.T[,1:nT]
				      }
				else  nT <- ncol(TT)



			} else {
				dbMEM.T=apply(as.matrix(COD.T),2,scale,center=TRUE,scale=TRUE)
				nT <- dim(dbMEM.T)[2]
				if(nT>=tt) stop("The number of temporal coding functions must be lower than t", call.=FALSE)
			  if(nrow(as.matrix(COD.T))!=dim(Y)[1]) stop("The number of rows in COD.T must be equal to the number of observations", call.=FALSE) 
			}
		}

#		if(s*tt-s-tt-nT*nS-1<=0 && model!="2" && model!="6a" && model!="6b") stop("Not enough degrees of freedom for testing interaction", call.=FALSE)

		if(print.res) {
			if(model!="2" && model!="3b" && model!="6b") { 
			cat(" Number of space coding functions =", nS,'\n\n')}            		
			if(model!="2" && model!="3a" && model!="6a") {
            		cat(" Number of time coding functions =", nT, "\n\n")}            		
		}
				

		# Generates space and time helmert contrasts
		A <- as.factor(rep(1:s,tt))
		B <- rep(1,s)
		for(i in 2:tt) B <- c(B,rep(i,s))
		B <- as.factor(B)
		HM <- model.matrix(~ A + B, contrasts = list(A="contr.helmert", B="contr.helmert"))
		HM.S <- as.matrix(HM[,2:s])
		HM.T <- as.matrix(HM[,(s+1):(s+tt-1)])
		
		test.STI = TRUE
		test.S = TRUE
		test.T = TRUE
		
		# Defines X (variables for the factor of interest) and W (covariables) for each test to be done
		if(model=="5") {
			XSTI <- dbMEM.S*dbMEM.T[,1]
			if(dim(dbMEM.T)[2]>1) for(j in 2:dim(dbMEM.T)[2]) XSTI <- cbind(XSTI,dbMEM.S*dbMEM.T[,j])
			XS <- HM.S
			XT <- HM.T
			if(print.res) {
				cat(" MODEL V: HELMERT CONTRAST FOR TESTING MAIN FACTORS. \n")
				cat("          SPACE AND TIME dbMEMs FOR TESTING INTERACTION.",'\n')
			}					
		} else if(model=="4") {
			XSTI <- dbMEM.S*dbMEM.T[,1]
			if(dim(dbMEM.T)[2]>1) for(j in 2:dim(dbMEM.T)[2]) XSTI = cbind(XSTI,dbMEM.S*dbMEM.T[,j])
			XS <- dbMEM.S
			XT <- dbMEM.T
			if(print.res) {
				cat(" MODEL IV: dbMEMs FOR BOTH SPACE AND TIME.",'\n')
			}		
		} else if(model=="3a") {
			XSTI <- dbMEM.S*HM.T[,1]
			if(tt>1) for(j in 2:(tt-1)) XSTI <- cbind(XSTI,dbMEM.S*HM.T[,j])
			XS <- dbMEM.S
			XT <- HM.T
			if(print.res) {
				cat(" MODEL IIIa: dbMEMs FOR SPACE AND HELMERT CONTRASTS FOR TIME.",'\n')
			}	
		} else if(model=="3b") {
			XSTI <- dbMEM.T*HM.S[,1]
			for(j in 2:(s-1)) XSTI <- cbind(XSTI,dbMEM.T*HM.S[,j])
			XS <- HM.S
			XT <- dbMEM.T
			if(print.res) {
				cat(" MODEL IIIb: HELMERT CONTRASTS FOR SPACE AND dbMEMs FOR TIME.",'\n')
			}		
		} else if(model=="7") {
			XSTI <- HM.S*HM.T[,1]
			if(dim(HM.T)[2]>1) for(j in 2:dim(HM.T)[2]) XSTI <- cbind(XSTI,HM.S*HM.T[,j])
			XS <- dbMEM.S
			XT <- dbMEM.T
			if(print.res) {
				cat(" MODEL VII: dbMEMs FOR BOTH SPACE AND TIME BUT HELMERT CONTRAST FOR INTERACTION.",'\n')
			}		
		} else if(model=="6a") {
			XS <- dbMEM.S
			for(j in 1:(tt-1)) XS <- cbind(XS,dbMEM.S*HM.T[,j])
			XT <- HM.T
			XSTI = NULL
			if(print.res) {
				cat(" MODEL VIa: NESTED MODEL.\n")
				cat("            TESTING FOR THE EXISTENCE OF SPATIAL STRUCTURE (COMMON OR SEPARATE)",'\n')
			}	
			test.STI = FALSE
		} else if(model=="6b") {
			XT <- dbMEM.T
			for(j in 1:(s-1)) XT <- cbind(XT,dbMEM.T*HM.S[,j])
			XS <- HM.S						
			XSTI = NULL
			if(print.res) {
				cat(" MODEL VIb: NESTED MODEL.\n")
				cat("            TESTING FOR THE EXISTENCE OF TEMPORAL STRUCTURE (COMMON OR SEPARATE).",'\n')
			}	
			test.STI = FALSE
		} else if(model=="2") {
			XS <- HM.S
			XT <- HM.T
			XSTI = NULL
			if(print.res) {
				cat(" MODEL II: HELMERT CONTRAST FOR SPACE AND TIME. NO INTERACTION TERM.",'\n')
			}	
			test.STI = FALSE
		} 		
			if(print.res) {
				cat("   Number of space variables =", dim(XS)[2],'\n')
				cat("   Number of time variables =", dim(XT)[2],'\n')
				if(test.STI) {
					cat("   Number of interaction variables =", dim(XSTI)[2],'\n')
					cat("   Number of residual degrees of freedom =", (s*tt-dim(XS)[2]-dim(XT)[2]-dim(XSTI)[2]-1),"\n\n")
					if((s*tt-dim(XS)[2]-dim(XT)[2]-dim(XSTI)[2]-1) <= 0) stop("Not enough residual degrees of freedom for testing. \nTry with a lower number of space or time variables.", call.=FALSE)
				} else {
					cat("   Number of residual degrees of freedom =", (s*tt-dim(XS)[2]-dim(XT)[2]-1),"\n\n")			
					if((s*tt-dim(XS)[2]-dim(XT)[2]-1) <= 0) stop("Not enough residual degrees of freedom for testing. \nTry with a lower number of space or time variables.", call.=FALSE)

				}
			}		
		
				
	res <- manovRDa(Y=Y,s=s,tt=tt,S.mat=XS,T.mat=XT,STI.mat=XSTI, Sfixed= Sfixed, Tfixed=Tfixed, S.test=test.S, T.test=test.T, STI.test=test.STI, model = model, nperm=nperm)
		
	if(test.STI==TRUE) {
			cat(' Interaction test:   R2 =', round(res$testSTI$R2, 4),'  F =',round(res$testSTI$F, 4),'  P(',nperm,'perm) =',res$testSTI$Prob,'\n')
		}
	if(test.S==TRUE) {
			cat(' Space test:         R2 =', round(res$testS$R2, 4),'  F =',round(res$testS$F, 4),'  P(',nperm,'perm) =',res$testS$Prob,'\n')
		}
	if(test.T==TRUE) {
			cat(' Time test:          R2 =', round(res$testT$R2, 4),'  F =',round(res$testT$F, 4),'  P(',nperm,'perm) =',res$testT$Prob,'\n')
		}
		
	})
	aa[3] <- sprintf("%2f",aa[3])
	if(print.res) {
		cat("-------------------------------------------------------\n")
		cat("      Time for this computation =",aa[3]," sec",'\n')
		cat("=======================================================\n\n")
	}
            class(res) <- "sti"
	invisible(res)
}




#' @rdname stimodels  
'quicksti' <- function(Y, S, Ti, nperm=999, alpha = 0.05, COD.S=NULL, COD.T=NULL,print.res=TRUE)
{
#	require(MASS)
#	require(adespatial)
	if(!is.logical(print.res)) {
		stop("Wrong operator; 'print.res' should be either 'FALSE' or 'TRUE'.", call.=FALSE)
	}

	aa <- system.time( {

	    # Sets the number and location of spatial and temporal points
		S <- as.matrix(S)
		if(dim(S)[1]==1 && dim(S)[2]==1) {
			s <- S[1,1]
			sitesX <- c(1:s)
		} else {
			s <- dim(S)[1]
			sitesX <- S
		}		
		Ti <- as.matrix(Ti)
		if(dim(Ti)[1]==1 && dim(Ti)[2]==1) {
			tt <- Ti[1,1]
			timesX <- c(1:tt) 
		} else {
			tt <- dim(Ti)[1]
			timesX <- Ti
		}

		
		# Total number of rows in matrix
		n <- s*tt			
		
		# Check response data file containing species data
		Y <- as.matrix(Y)
		p <- dim(Y)[2]
		if(dim(Y)[1] != n) stop("The number of rows in species file is not (S x Ti).", call.=FALSE) 

		# Center response data
		Y <- scale(Y, center=TRUE, scale=FALSE)
		
				
		if(print.res) {
			cat("=========================================================\n")
			cat("        Space-time ANOVA without replicates\n")
			cat("                                                  \n")
			cat("  Pierre Legendre, Miquel De Caceres, Daniel Borcard\n")
			cat("---------------------------------------------------------\n\n")

			cat(" Number of space points (s) =", s,'\n')
			cat(" Number of time points (tt) =", tt,'\n')
			cat(" Number of observations (n = s*tt) =", n,'\n')
			cat(" Number of response variables (p) =", p,'\n')
			cat(" Significance level for the interaction test (alpha) =", alpha,'\n','\n')
			
		}

		# Generates spatial dbMEMs if not given by user
		if(is.null(COD.S)) {			
			if(print.res) cat(" Computing dbMEMs to code for space\n")
			dbMEM.S.tmp <- dbmem(sitesX, MEM.autocor="positive")
			SS <- as.matrix(dbMEM.S.tmp)
			   nS <- ncol(SS)
			dbMEM.S.thresh <- give.thresh(dist(sitesX))
			if(print.res) cat(" Truncation level for space dbMEMs =", dbMEM.S.thresh, "\n")
			dbMEM.S <- SS
			for(j in 2:tt) dbMEM.S = rbind(dbMEM.S,SS)
		} else {
			dbMEM.S <- apply(as.matrix(COD.S),2,scale,center=TRUE,scale=TRUE)
			nS <- dim(dbMEM.S)[2]
			if(nS >= s) stop("The number of spatial coding functions must be lower than S.", call.=FALSE) 
			if(nrow(as.matrix(COD.S))!=dim(Y)[1]) stop("The number of rows in COD.S must be equal to the number of observations.", call.=FALSE) 
		}
			
		# Generate dbMEM variables for time if not given by user
		if(is.null(COD.T)) {
			nT <- trunc(tt/2)
			if(print.res) cat(" Computing dbMEMs to code for time\n")
			dbMEM.T.tmp <- dbmem(timesX, MEM.autocor="positive")
			TT <- as.matrix(dbMEM.T.tmp)
			   nT <- ncol(TT)
			dbMEM.T.thresh <- give.thresh(dist(timesX))
			if(print.res) cat(" Truncation level for time dbMEMs =", dbMEM.T.thresh, "\n\n")
			T.temp <- TT[1,]
			for(i in 2:s) T.temp <- rbind(T.temp,TT[1,])
			dbMEM.T <- as.matrix(T.temp)
			for(j in 2:tt) {
				T.temp <- TT[j,]
				for(i in 2:s) T.temp <- rbind(T.temp,TT[j,])
				dbMEM.T <- as.matrix(rbind(dbMEM.T,T.temp))
			}
		} else {
			dbMEM.T <- apply(as.matrix(COD.T),2,scale,center=TRUE,scale=TRUE)
			nT <- dim(dbMEM.T)[2]
			if(nT >= tt) stop("The number of temporal coding functions must be lower than Ti.", call.=FALSE)
			if(nrow(as.matrix(COD.T))!=dim(Y)[1]) stop("The number of rows in COD.T must be equal to the number of observations.", call.=FALSE) 
		}
		
		if(s*tt-s-tt-nT*nS-1<=0) stop("Not enough degrees of freedom for testing interaction.", call.=FALSE)


		if(print.res) {
			cat(" Number of space coding functions =", nS,'\n')
            		cat(" Number of time coding functions =", nT, "\n\n")
		}

		# Generates space and time helmert contrasts
		A <- as.factor(rep(1:s,tt))
		B <- rep(1,s)
		for(i in 2:tt) B <- c(B,rep(i,s))
		B <- as.factor(B)
		HM <- model.matrix(~ A + B, contrasts = list(A="contr.helmert", B="contr.helmert"))
		HM.S <- as.matrix(HM[,2:s])
		HM.T <- as.matrix(HM[,(s+1):(s+tt-1)])

		# Test significance for interaction effect 
		#
		# Defines X (variables for the factor of interest) and W (covariables) for the space-time test
		#
		XSTI <- dbMEM.S*dbMEM.T[,1]
		if(dim(dbMEM.T)[2]>1) for(j in 2:dim(dbMEM.T)[2]) XSTI <- cbind(XSTI,dbMEM.S*dbMEM.T[,j])
		if(print.res) {
			cat("------------------------------------------\n")
			cat(" Testing space-time interaction (model 5)\n")
			cat("------------------------------------------\n\n")
			cat("   Number of space variables =", dim(HM.S)[2],'\n')
			cat("   Number of time variables =", dim(HM.T)[2],'\n')
			cat("   Number of interaction variables =", dim(XSTI)[2],'\n')
			cat("   Number of residual degrees of freedom =", (s*tt-dim(XSTI)[2]-dim(HM.S)[2]-dim(HM.T)[2]-1),"\n\n")
			if((s*tt-dim(XSTI)[2]-dim(HM.S)[2]-dim(HM.T)[2]-1) <= 0) stop("Not enough residual degrees of freedom for testing. Try with a lower number of space or time variables.", call.=FALSE)
		}
				
		res <- manovRDa(Y=Y,s=s,tt=tt,S.mat=HM.S,T.mat=HM.T,STI.mat=XSTI, S.test=TRUE, T.test=TRUE, STI.test=TRUE, model = "5", nperm=nperm)
		
		if(print.res) {
			cat(' Interaction test:  R2 =', round(res$testSTI$R2, 4),'  F =',round(res$testSTI$F, 4),'  P(',nperm,'perm) =',res$testSTI$Prob,"\n\n")
		}
						
		if(res$testSTI$Prob<=alpha) {
			XS <- dbMEM.S
			for(j in 1:(tt-1)) XS = cbind(XS,dbMEM.S*HM.T[,j])
			XT <- HM.T
			if(print.res) {
				cat("----------------------------------------------------------------------\n")
				cat(" Testing for the existence of separate spatial structures (model 6a)\n")
				cat("----------------------------------------------------------------------\n\n")
				cat("   Number of space variables =", dim(XS)[2],'\n')
				cat("   Number of time variables =", dim(XT)[2],'\n')
				cat("   Number of residual degrees of freedom =", (s*tt-dim(XS)[2]-dim(XT)[2]-1),"\n\n")
			}			
			if((s*tt-dim(XS)[2]-dim(XT)[2]-1)==0) {
				cat("   Not enough degrees of freedom to perform test with model 6b.\n   Try with a lower number of coding variables for space.\n")
			} else {
				res2 <- manovRDa(Y=Y,s=s,tt=tt,S.mat=XS,T.mat=XT,STI.mat=NULL, S.test=TRUE, T.test=FALSE, STI.test=FALSE, model="6a", nperm=nperm)
				res$testS <- res2$testS
				cat(' Space test:  R2 =', round(res$testS$R2,4),'  F =',round(res$testS$F,4),'  P(',nperm,'perm) =',res$testS$Prob,"\n\n")
			}
			XT <- dbMEM.T
			for(j in 1:(s-1)) XT <- cbind(XT,dbMEM.T*HM.S[,j])
			XS <- HM.S
			if(print.res) {
				cat("-----------------------------------------------------\n")
				cat(" Testing for separate temporal structures (model 6b)\n")
				cat("-----------------------------------------------------\n\n")
				cat("   Number of space variables =", dim(XS)[2],'\n')
				cat("   Number of time variables =", dim(XT)[2],'\n')
				cat("   Number of residual degrees of freedom =", (s*tt-dim(XS)[2]-dim(XT)[2]-1),"\n\n")
			}	
			if((s*tt-dim(XS)[2]-dim(XT)[2]-1)==0) {
				cat("   Not enough degrees of freedom to perform test with model 6b.\n   Try with a lower number of coding variables for time.\n")
			} else {
				res3 <- manovRDa(Y=Y,s=s,tt=tt,S.mat=XS,T.mat=XT,STI.mat=NULL, S.test=FALSE, T.test=TRUE, STI.test=FALSE, model="6b", nperm=nperm)
				res$testT <- res3$testT
				cat(' Time test:   R2 =', round(res$testT$R2,4),'  F =',round(res$testT$F,4),'  P(',nperm,'perm) =',res$testT$Prob,"\n\n")			}		
		} else {
			if(print.res) {
				cat("---------------------------------------------------------------------\n")
				cat(" Testing for common spatial and common temporal structures (model 5)\n")
				cat("---------------------------------------------------------------------\n\n")
				cat(' Space test:   R2 =', round(res$testS$R2,4),'  F =',round(res$testS$F,4),'  P(',nperm,'perm) =',res$testS$Prob,'\n')
				cat(' Time test:    R2 =', round(res$testT$R2,4),'  F =',round(res$testT$F,4),'  P(',nperm,'perm) =',res$testT$Prob,'\n')
			}	
		} 		
		
	})
	
	aa[3] <- sprintf("%2f",aa[3])
	if(print.res) {
		cat("---------------------------------------------------------\n")
		cat("      Time for this computation =",aa[3]," sec",'\n')
		cat("=========================================================\n\n")
	}
	invisible(res)
}




'manovRDa' <- function(Y, s, tt, S.mat=NULL, T.mat=NULL, STI.mat=NULL, Sfixed=TRUE, Tfixed=TRUE, S.test=TRUE, T.test=TRUE, STI.test=TRUE, model = "5", nperm=999)
#
# ARGUMENTS:
#
# Y     		site-by-species data table (assumes row blocks corresponding to times, i.e. within each block 
#				all sites occur)
# s     		number of spatial points
# tt     		number of time campaigns
# S.mat 		a matrix of spatial variables.
# T.mat 		a matrix of temporal variables.
# STI.mat 		a matrix of interaction variables.
# Sfixed		logical: is Factor space fixed, or not (if FALSE, it is considered a random factor)
# Tfixed		logical: is Factor time fixed, or not (if FALSE, it is considered a random factor)
#
# S.test		logical: should space be tested (TRUE) or not (FALSE)?
# T.test		logical: should time be tested (TRUE) or not (FALSE)?
# STI.test		logical:  should interaction be tested (TRUE) or not (FALSE)?
#
# model 		space-time model (used to identify nested models "6a" and "6b")
#
# nperm 		number of permutations to be done
#
#===============================================================================
#
{
# require(MASS)

n <- nrow(Y)
p <- ncol(Y)
a <- ncol(S.mat)
b <- ncol(T.mat)

if(!is.null(STI.mat)) {cc <- ncol(STI.mat)}
else {cc = 0}

A <- S.mat
B <- T.mat
if(!is.null(STI.mat)) AxB <- STI.mat

# Compute projector of A and Yfit.A
invA <- ginv(t(A) %*% A)
projA <- A %*% invA %*% t(A)
Yfit.A <- projA %*% Y

# Compute projector of B and Yfit.B
invB <- ginv(t(B) %*% B)
projB <- B %*% invB %*% t(B)
Yfit.B <- projB %*% Y

# Compute projector of AxB and Yfit.AxB
if(!is.null(STI.mat)) {
	invAxB <- ginv(t(AxB) %*% AxB)
	projAxB <- AxB %*% invAxB %*% t(AxB)
	Yfit.AxB <- projAxB %*% Y
}

# Create a "compound matrix" to obtain R-square and adjusted R-square
if(!is.null(STI.mat)) {
	ABAxB <- cbind(A,B,AxB)
} else {
	ABAxB <- cbind(A,B)
}

# Compute projector of ABAxB and Yfit.ABAxB
invABAxB <- ginv(t(ABAxB) %*% ABAxB)
projABAxB <- ABAxB %*% invABAxB %*% t(ABAxB)
Yfit.ABAxB <- projABAxB %*% Y

# Compute Sums of squares (SS) and Mean squares (MS)
SS.Y <- sum(Y*Y)
SS.Yfit.ABAxB <- sum(Yfit.ABAxB*Yfit.ABAxB)
SS.Yfit.A <- sum(Yfit.A*Yfit.A)
SS.Yfit.B <- sum(Yfit.B*Yfit.B)
if(!is.null(STI.mat)) SS.Yfit.AxB <- sum(Yfit.AxB*Yfit.AxB)
MS.A <- SS.Yfit.A/a
MS.B <- SS.Yfit.B/b
if(!is.null(STI.mat)) MS.AxB <- SS.Yfit.AxB/cc
MS.Res <- (SS.Y-SS.Yfit.ABAxB)/(n-(a+b+cc)-1)


if(STI.test==TRUE) { # Test interaction (unrestricted permutations)
	nPGE.AxB = 1
	Fref.AxB <- MS.AxB/MS.Res
	vec <- c(1:n)
	for(i in 1:nperm) {
		YPerm <- Y[restrictedPerm(nobs.block=s, nblock=tt, n, restPerm=0, vec),]
		YhatPerm.AxB <- projAxB %*% YPerm
		SS.YhatPerm.AxB <- sum(YhatPerm.AxB*YhatPerm.AxB)
		MS.Perm.AxB <- SS.YhatPerm.AxB/cc
	
		YhatPerm.ABAxB <- projABAxB %*% YPerm
		SS.YhatPerm.ABAxB <- sum(YhatPerm.ABAxB*YhatPerm.ABAxB)

		MS.Perm.Res <- (SS.Y-SS.YhatPerm.ABAxB)/(n-(a+b+cc)-1)

		Fper.AxB <- MS.Perm.AxB/MS.Perm.Res
	 	if(Fper.AxB >= Fref.AxB) nPGE.AxB <- nPGE.AxB+1
	}
	P.AxB <- nPGE.AxB/(nperm+1)

	R2 <- SS.Yfit.AxB/SS.Y
	R2a <- 1-((n-1)/(n-dim(STI.mat)[2]-1))*(1-R2)
	
	testSTI <- list(MS.num=MS.AxB, MS.den=MS.Res, R2=R2, R2.adj=R2a, F=Fref.AxB, Prob=P.AxB)
} else {
	testSTI <- NULL
}


if(S.test==TRUE) { # Test factor A (space) using restricted permutations within time blocks
	nPGE.A=1
	if(Tfixed==FALSE && !is.null(STI.mat)) { # Time random factor in crossed design with interaction
		Fref.A <- MS.A/MS.AxB
		MS.den = MS.AxB
	} else if(Tfixed==FALSE && model=="6b") { # Time random factor in nested design
		Fref.A <- MS.A/MS.B
		MS.den <- MS.B
	} else {
		Fref.A=MS.A/MS.Res
		MS.den = MS.Res
	}
	vec <- c(1:n)
	for(i in 1:nperm) {
		YPerm <- Y[restrictedPerm(nobs.block=s,nblock=tt, n,restPerm=1, vec),]
		YhatPerm.A <- projA %*% YPerm
		SS.YhatPerm.A <- sum(YhatPerm.A*YhatPerm.A)
		MS.Perm.A <- SS.YhatPerm.A/a

		if(Tfixed==FALSE && !is.null(STI.mat)) { # Time random factor in crossed design with interaction
			YhatPerm.AxB <- projAxB %*% YPerm
			SS.YhatPerm.AxB <- sum(YhatPerm.AxB*YhatPerm.AxB)
			MS.Perm.AxB <- SS.YhatPerm.AxB/cc
			Fper.A <- MS.Perm.A/MS.Perm.AxB
		} else if(Tfixed==FALSE && model=="6b") { # Time random factor in nested design
			YhatPerm.B <- projB %*% YPerm
			SS.YhatPerm.B <- sum(YhatPerm.B*YhatPerm.B)
			MS.Perm.B <- SS.YhatPerm.B/b
			Fper.A <- MS.Perm.A/MS.Perm.B
		} else {
			YhatPerm.ABAxB <- projABAxB %*% YPerm
			SS.YhatPerm.ABAxB <- sum(YhatPerm.ABAxB*YhatPerm.ABAxB)
			MS.Perm.Res <- (SS.Y-SS.YhatPerm.ABAxB)/(n-(a+b+cc)-1)
			Fper.A <- MS.Perm.A/MS.Perm.Res
		}
		if(Fper.A >= Fref.A) nPGE.A <- nPGE.A+1
	}
	P.A <- nPGE.A/(nperm+1)

	R2 <- SS.Yfit.A/SS.Y
	R2a <- 1-((n-1)/(n-a-1))*(1-R2)
		
	testS <- list(MS.num=MS.A, MS.den=MS.den, R2=R2, R2.adj=R2a, F=Fref.A, Prob=P.A)
} else {
	testS <- NULL
}

if(T.test==TRUE) { # Test factor B (time) using restricted permutations within time blocks
	nPGE.B = 1
	if(Sfixed==FALSE && !is.null(STI.mat)) { # Space random factor in crossed design with interaction
		Fref.B <- MS.B/MS.AxB
		MS.den <- MS.AxB
	} else if(Sfixed==FALSE && model=="6a") { # Space random factor in nested design
		Fref.B <- MS.B/MS.A
		MS.den <- MS.A
	} else {
		Fref.B <- MS.B/MS.Res
		MS.den <- MS.Res
	}

	vec <- c(1:n)
	for(i in 1:nperm) {
		YPerm <- Y[restrictedPerm(nobs.block=s,nblock=tt, n,restPerm=2,vec),]
		YhatPerm.B <- projB %*% YPerm
		SS.YhatPerm.B <- sum(YhatPerm.B*YhatPerm.B)
		MS.Perm.B <- SS.YhatPerm.B/b

		if(Sfixed==FALSE && !is.null(STI.mat)) { # Space random factor in crossed design with interaction
			YhatPerm.AxB <- projAxB %*% YPerm
			SS.YhatPerm.AxB <- sum(YhatPerm.AxB*YhatPerm.AxB)
			MS.Perm.AxB <- SS.YhatPerm.AxB/cc
			Fper.B <- MS.Perm.B/MS.Perm.AxB
		} else if(Sfixed==FALSE && model=="6a") { # Space random factor in nested design
			YhatPerm.A <- projA %*% YPerm
			SS.YhatPerm.A <- sum(YhatPerm.A*YhatPerm.A)
			MS.Perm.A <- SS.YhatPerm.A/a
			Fper.B <- MS.Perm.B/MS.Perm.A
		} else {
			YhatPerm.ABAxB <- projABAxB %*% YPerm
			SS.YhatPerm.ABAxB <- sum(YhatPerm.ABAxB*YhatPerm.ABAxB)
			MS.Perm.Res <- (SS.Y-SS.YhatPerm.ABAxB)/(n-(a+b+cc)-1)
			Fper.B <- MS.Perm.B/MS.Perm.Res
		}	
		if(Fper.B >= Fref.B) nPGE.B <- nPGE.B+1
	}
	P.B <- nPGE.B/(nperm+1)

	R2 <- SS.Yfit.B/SS.Y
	R2a <- 1-((n-1)/(n-b-1))*(1-R2)
	
	testT <- list(MS.num=MS.B, MS.den=MS.den, R2=R2, R2.adj=R2a, F=Fref.B, Prob=P.B)
} else {
	testT <- NULL
}

return(list(testSTI = testSTI, testS = testS, testT=testT))
}

'restrictedPerm' <- function(nobs.block, nblock, n, restPerm, vec)
#
# restPerm == 0: Unrestricted permutation. 
#
# restPerm == 1: Restricted permutation of the observations within each block.
# The data are arranged by blocks, with all observations forming a block 
# placed in adjacent positions in the file. Example:
# BLOCK-1: obs1, obs2, ... obs-nobs.block; BLOCK-2: obs1, obs2, ... obs-nobs.block; etc.
#
# restPerm == 2: Restricted permutation of the observations across blocks (within the blocks formed by each nth 
# observation).
#
# Vector 'vec' contains the initial order of the objects, e.g. vec=c(1:n).
# At the end of the function, it gives the permuted order of the objects.
#
# Examples:  toto0 <- restrictedPerm(6,4,24,0,c(1:24))
#            toto1 <- restrictedPerm(6,4,24,1,c(1:24))
#
#                                       Pierre Legendre, January 2006
#			 Miquel De Caceres, February 2009
{
	if(restPerm == 0) { 
		vec <- sample(vec[1:n],n)
	} else if(restPerm==1) {
		for(j in 1:nblock) {
			i1 <- nobs.block*(j-1)+1
			i2 <- nobs.block*j
			vec[i1:i2] <- sample(vec[i1:i2],nobs.block)
		}
	} else {
		vecT <- vec
		for(j in 1:nobs.block) {
			ind <- sample(1:nblock)
			for(i in 1:nblock) {
				vec[nobs.block*(i-1)+j] <- vecT[nobs.block*(ind[i]-1)+j]
			}
		}		
	}

	return(vec)
}
