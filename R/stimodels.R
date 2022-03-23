#' Space-time interaction in ANOVA without replication
#' 
#' Function \code{stimodels} performs two-way ANOVA to test space-time
#' interaction (STI) without replicates using one among a set of possible models 
#' described in Legendre et al. (2010).
#' Function \code{quicksti} allows performing space-time ANOVA in a simplified
#' way. In many models, degrees of freedom are saved by coding space and/or
#' time parsimoniously using distance-based Moran Eigenvector Maps (dbMEM; 
#' Borcard & Legendre 2002; Dray et al. 2006).
#' 
#' @details The 'stimodels' and 'quicksti' functions should only be used 
#' (1) when each site has been sampled during each survey, with no missing data, 
#' and (2) when there are no replicate observations of the space-time points. When 
#' there is replication, use a regular permutational Manova function such as adonis2.
#' 
#' When the sites form a one-dimensional spatial transect, or a meandering line 
#' such as the course of a river, with regular sampling intervals, and the time series 
#' has fairly equal time intervals, one can use the S and Ti arguments to indicate 
#' the number of sites in space and the number of surveys along time. The order of 
#' the sites in each temporal block of the input data file will be taken to correspond 
#' to the spatial positions of the sites along the transect (from 1 to S), and the 
#' order of the time blocks in the data file will be taken to indicate the temporal 
#' order of the surveys (from 1 to Ti). The function will then compute dbMEM 
#' eigenfunctions corresponding to the spatial and temporal positions of the data 
#' rows in the input data file as if they were on straight lines.
#' 
#' When the sites do not form a regular, one-dimensional spatial transect, one must 
#' provide a file of spatial coordinates of the sites to argument S. Similarly, when 
#' the time series has unequal time intervals, one must provide a file of temporal 
#' coordinates of the surveys to argument Ti.
#' 
#' Alternatively, one can use arguments COD.S or COD.T to provide a matrix of Helmert 
#' contrasts to be used in the analysis in place of dbMEM eigenfunctions. One can do 
#' that, for example, when there are only a few surveys along time and it would not 
#' be useful to represent these few surveys by dbMEM eigenfunctions. That matrix can 
#' have the class "matrix" or "numeric"; they both work in functions stimodels and 
#' quicksti. Arguments COD.S and COD.T can also be used to provide matrices 
#' containing other types of eigenfunctions, for example AEM eigenfunctions, to be 
#' used instead of dbMEM matrices computed by stimodels or quicksti. However, do not 
#' code both the space and time factors as Helmert contrasts; that would not leave 
#' residual degrees of freedom for the test of the interaction.
#' 
#' In \code{stimodels}, tests for space-time interaction and space or time main
#' effects are conducted using one of the following models: 
#' \itemize{
#' \item Model 2 - Space and Time are coded using Helmert contrasts for the main
#' effects. The interaction cannot be tested. 
#' \item Model 3a - Space is coded using dbMEM variables whereas Time is coded 
#' using Helmert contrasts. 
#' \item Model 3b - Space is coded using Helmert contrasts whereas Time is coded
#' using dbMEM variables.
#' \item Model 4 - Both Space and Time are coded using dbMEM variables; the 
#' interaction is coded as the Hadamard (or elementwise) product of the 
#' space-coding by the time-coding dbMEMs.
#' \item Model 5 - Space and Time are coded using Helmert contrasts for the main
#' factor effects; the interaction is coded as the Hadamard product of the 
#' space-coding by the time-coding dbMEM variables. 
#' \item Model 6a - Nested model. Testing for the existence of spatial
#' structure (common or separate) using dbMEM (or other) variables coding for Space.
#' \item Model 6b - Nested model. Testing for the existence of temporal structure
#' (common or separate) using dbMEM (or other) variables coding for Time. 
#' \item Model 7 - Space and Time are coded using dbMEM variables for the main
#' factor effects, but they are coded using Helmert contrasts for the interaction 
#' term (not recommended). }
#' 
#' With Models 2, 6a and 6b, the interaction test is not available. 
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
#' @param Y Site-by-species response data table. Assumes row blocks
#' corresponding to times, i.e. within each block all sites are provided, 
#' always in the same order.
#' @param S Number of spatial points (when they are aligned on a transect or a
#' time series and equispaced) or a matrix of spatial coordinates (when the
#' sites are on a two-dimensional surface or on a line but very irregularly
#' spaced).
#' @param Ti Number of time campaigns (when equispaced) or a matrix (a
#' vector) of temporal coordinates (when the time campaigns are very 
#' irregularly spaced).
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
#' number of columns must be lower than \code{S} and the number of rows must
#' be equal to the number of rows in \code{Y}. – Do not use full coding of S 
#' (binary variables or Helmert contrasts) in COD.S to code for the factor to 
#' be tested, space, in Model 6a: there would be no residual degrees of freedom 
#' left for the test.
#' @param COD.T Temporal coding functions to be used instead of dbMEM. The
#' number of columns must be lower than \code{Ti} and the number of rows must
#' be equal to the number of rows in \code{Y}. – Do not use full coding of Ti 
#' (binary variables or Helmert contrasts) in COD.T to code for the factor to
#' be tested, time, in Model 6b: there would be no residual degrees of freedom 
#' left for the test.
#' @param save.X If TRUE, the explanatory bloc-diagonal matrix constructed for
#' model 6a or 6b is saved in the output list with the label X.
#' @param print.res If TRUE, displays the results and additional information
#' onscreen (recommended). If FALSE, only R2, F and P are printed onscreen.
#' 
#' @return A list containing the following results: \itemize{ 
#' \item{testS}{ An object with the result of the space effect test, including
#' the mean squares for the F numerator (\code{MS.num}), the mean squares for
#' the F denominator (\code{MS.den}), the proportion of explained variance
#' (\code{R2}), the adjusted proportion of explained variance (\code{R2.adj}),
#' the F statistics (\code{F}) and its p-value computed from a permutation test
#' (\code{Prob}). } 
#' \item{testT}{ An object with the result of the time effect test, like \code{testS}.} 
#' \item{teststi}{ An object with the result of the space-time interaction test, 
#' like \code{testS}.} 
#' \item{X.matrix}{ The bloc-diagonal explanatory matrix used in test of model 6a or 6b}
#' }
#'
#' @author Pierre Legendre \email{pierre.legendre@@umontreal.ca}, Miquel De Caceres 
#' and Daniel Borcard
#' @seealso \code{\link{trichoptera}}
#'
#' @references 
#' Borcard, D. & P. Legendre. 2002. All-scale spatial analysis of ecological data by 
#' means of principal coordinates of neighbour matrices. Ecological Modelling 153: 51-68.
#' https://doi.org/10.1016/S0304-3800(01)00501-4.
#'
#' Dray, S., P. Legendre & P. R. Peres-Neto. 2006. Spatial modelling: a
#' comprehensive framework for principal coordinate analysis of neighbour
#' matrices (PCNM). Ecological Modelling 196: 483-493. 
#' https://doi.org/10.1016/j.ecolmodel.2006.02.015.
#'
#' Legendre, P. & D. Borcard. 2018. Box-Cox-chord transformations for community 
#' composition data prior to beta diversity analysis. Ecography 41: 1820–1824. 
#' https://doi.org/10.1111/ecog.03498.
#' 
#' Legendre, P., M. De Caceres & D. Borcard. 2010. Community surveys through space 
#' and time: testing the space-time interaction in the absence of replication. 
#' Ecology 91: 262-272. https://doi.org/10.1890/09-0199.1.
#'
#' @keywords multivariate spatial
#' 
#' @examples
#' # The "trichoptera" data set contains the abundances of 56 Trichoptera species captured 
#' # during 100 consecutive days in 22 emergence traps along a stream. The 22 traps 
#' # (sites) form a regular transect, with geographic positions 1 to 22. The original 
#' # daily data collected at each site were pooled into 10 survey periods for the study 
#' # of Legendre et al. (2010) in order to reduce the very high proportion of zeros in the 
#' # original data matrix. Order of the observations in the data set: the 22 traps (sites) 
#' # are nested within the survey weeks, as required by the 'stimodels' and 'quicksti' 
#' # functions.
#' 
#' data(trichoptera)
#' 
#' # log-transform the species data, excluding the Site and Date colums
#' 
#' trich.log <- log1p(trichoptera[,-c(1,2)]) 
#' 
#' # A log-chord transformation (Legendre & Borcard 2018) would also be appropriate for 
#' # these data: trich.tr <- decostand(log1p(trichoptera[,-c(1,2)]), method="norm")
#' 
#' # Example 1. Compute the space-time interaction test using model 5. By specifying the  
#' # number of sites (traps), the sofware assumes that they form a regular transect with  
#' # equispaced sites. Important note to users – In real analyses, use more than 99 
#' # permutations.
#' 
#' out.1 <- stimodels(trich.log, S=22, Ti=10, model="5", nperm=99)
#' 
#' # The interaction is significant. Because of that, test results for the main effects, 
#' # space and time, obtained with model 5, cannot be interpreted. Tests of spatial 
#' # variation can be done for individual times using simple RDA against dbMEM.S  
#' # variables. Likewise, tests of temporal variation can be done for individual sites  
#' # using simple RDA against dbMEM.T variables. A global test of the hypothesis that none 
#' # of the times shows a significant spatial structure can be done with model 6a. For a  
#' # global test of temporal structure at the various sites, use model 6b.
#' 
#' \dontrun{     # Code not run during CRAN software tests
#'
#' # Example 2. Run space-time analysis with global tests for main effects after testing
#' # the interaction, which is significant in this example 
#' 
#' out.2 <- quicksti(trich.log, S=22, Ti=10, nperm=999)
#'
#' # Since the interaction is significant, function 'quicksti' will carry out the 
#' # tests of existence of a spatial (at least in one of the time periods) and temporal 
#' # (at least at one of the sites) structures using models 6a and 6b, respectively.
#' 
#' # 3. Run space-time analysis for two time blocks only, i.e. time periods 6 and 7,   
#' # then time periods 8 and 9.
#' 
#' # Example 3.1. Time periods 6 and 7. The interaction is not significant. In that case,  
#' # function 'quicksti' carries out the tests of the main effects using model 5.
#' 
#' out.3 <- quicksti(trich.log[111:154,], S=22, Ti=2, nperm=999)
#'
#' # Example 3.2. Time periods 8 and 9. The interaction is significant. In that case,  
#' # 'quicksti' carries out the tests of the spatial effects using model 6a. Model 6b 
#' # cannot proceed with the test of the temporal effect because Ti=2. An explanation is 
#' # printed in the output list.
#' 
#' out.4 <- quicksti(trich.log[155:198,], S=22, Ti=2, nperm=999)
#' 
#' # 4. Illustrations of the use of 'COD.S' and 'COD.T' in STI analysis
#'
#' # The following examples illustrate how to use other representations of the spatial or 
#' # temporal relationships among observations, through arguments 'COD.S' and 
#' # 'COD.T' of functions 'stimodels' and 'quicksti'. The trichoptera data 
#' # are used again.
#'
#' # Example 4.1. Explicitly compute dbMEMs for the spatial structure along the regular 
#' # transect, using function 'dbmem' of adespatial, and provide it to 'stimodels' 
#' # or 'quicksti' as argument 'COD.S'. The dbMEMs must first be computed on the
#' # transect, then repeated (Ti-1) times to provide Ti repeats in total.
#'
#' dbMEM.S1 <- as.matrix(dbmem(1:22))
#' dbMEM.S10 <- dbMEM.S1
#' for(j in 2:10) dbMEM.S10 <- rbind(dbMEM.S10, dbMEM.S1)
#' out.5 <- stimodels(trich.log, S=22, Ti=10, model="5", COD.S=dbMEM.S10, nperm=999)
#'
#' # Results should be identical to those in output file out.1 of Example 1, except for 
#' # P-values which can vary slightly.
#'
#' # Example 4.2. Assume now that the sampling sites have irregular positions, as 
#' # described by the following matrix of geographic coordinates 'xy.trich'. Provide 
#' # this matrix to argument S of 'stimodels'
#'
#' xy.trich = matrix(c(1:5,11:15,21:25,31:35,41,42,rep(c(1,2),11)),22,2)
#' plot(xy.trich, asp=1)   # Plot a quick map of the site positions
#' out.6 <- stimodels(trich.log, S=xy.trich, Ti=10, model="5", nperm=999)
#'
#' # Example 4.3. Compute a matrix of dbMEMs for the sites. The coding matrix provided to 
#' # argument 'COD.S' must contain repeated dbMEM.S codes because that matrix must have 
#' # the same number of rows as matrix Y. Construct coding matrix dbMEM.S10 containing the 
#' # dbMEM.S codes repeated 10 times. 
#'
#' dbMEM.S1 <- as.matrix(dbmem(xy.trich))
#' dbMEM.S10 = dbMEM.S1
#' for(i in 1:9) dbMEM.S10 <- rbind(dbMEM.S10, dbMEM.S1)
#' out.7 <- stimodels(trich.log, S=22, Ti=10, model="5", COD.S=dbMEM.S10, nperm=999)
#'
#' # Compare the results with those obtained in the output file out6, example 4.2.
#'
#' # Careful: If an analysis requires a dbMEM coding matrix for 'COD.T', the dbMEM.T    
#' # codes must follow the required data arrangement: sites must be nested within times.
#' # The following function can be used to construct a dbMEM.T matrix.
#' 
#' MEM.T <- function(s, tt, coord.T=NULL)
#'   # Documentation of function MEM.T –
#'	 # Generate a matrix of temporal eigenfunctions for input into stimodels, 
#'	 # with sites nested within times.
#'	 # Arguments –
#'	 # s : number of space points (sites)
#'	 # tt : number of time points
#'	 # coord.T : optional matrix or vector giving the time point coordinates
#'   {
#'	  n <- s*tt
#'	  if(is.null(coord.T)) coord.T <- as.matrix(1:tt)
#'	  MEM.TT <- as.matrix(dbmem(coord.T))
#'	  dbMEM.T <- matrix(0,n,ncol(MEM.TT))    # Empty matrix to house dbMEM.T		
#'	  beg.x <- seq(1, n, by=s)
#'	  for(i in 1:tt) { # Fill tt blocks of rows with identical MEM.TT values
#'		  for(j in 1:s) dbMEM.T[(beg.x[i]+j-1),] <- MEM.TT[i,]
#'		  }
#'	  dbMEM.T
#'   }
#' 
#' # Example of use of function MEM.T
#'  
#' dbMEM.T <- MEM.T(s=6, tt=5)
#' # Check the size of the dbMEM.T output matrix
#' dim(dbMEM.T)
#' 
#' }   # End of code not run during CRAN software tests
#' 
#' @importFrom MASS ginv
#' @importFrom stats model.matrix
#' @aliases stimodels quicksti
#' @export stimodels quicksti
#' @rdname stimodels  
#'   

'stimodels' <- function(Y, S, Ti,  model="5", nperm=999, nS=-1, nT=-1, Sfixed=TRUE, 
                        Tfixed=TRUE, COD.S=NULL, COD.T=NULL, save.X=FALSE, print.res=TRUE)
{
    if(!is.logical(print.res)) {
        stop("Wrong operator; 'print.res' should be either 'FALSE' or 'TRUE'",
             call.=FALSE)
    }
    if(model!="2" && model!="3a" && model!="3b" && model!="4" && model!="5" 
       && model!="6a" && model!="6b" && model!="7") {
        stop(paste("Unrecognized model ",model,
                   "; 'model' should be '2', '3a', '3b', '4', '5', '6a','6b' or '7'.", 
                   sep=""), call.=FALSE)
    }
    # Checking that the data represent abundances or abundance-like data
    if (any(Y < 0))
        stop("Data contain negative values\n", call.=FALSE)
    if (any(is.na(Y)))
        stop("Data contain 'NA' values\n", call.=FALSE)
    
    skip.Manova <- FALSE
    res <- NA
    
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
        
        if(dim(Y)[1] != n) stop("The number of rows in species file is not (S x Ti)",
                                call.=FALSE) 
        
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
        
        # Generates space and time helmert contrasts
        A <- as.factor(rep(1:s,tt))
        B <- rep(1,s)
        for(i in 2:tt) B <- c(B,rep(i,s))
        B <- as.factor(B)
        HM <- model.matrix(~ A+B, contrasts = list(A="contr.helmert", B="contr.helmert"))
        HM.S <- as.matrix(HM[,2:s])            # Helmert contrasts for Space, no intercept
        HM.T <- as.matrix(HM[,(s+1):(s+tt-1)]) # Helmert contrasts for Time, no intercept
        
        # Generate dbMEM variables for Space (if necessary)
        if(model=="3a"|| model=="6a" || model=="7"|| model=="4"|| model=="5") {
            if(is.null(COD.S)) {	# Generate spatial dbMEMs if not given by user
                if(print.res) cat(" Computing dbMEMs to code for space\n")
                if(s==2) {
                    dbMEM.S <- as.matrix(rep(c(-1,1),tt))
                    #
                    if(model=="6a") { 
                        cat("Model 6a requires that 's' be larger than 2. When s=2, full",
                            "coding of the site locations by a binary variable or Helmert",
                            "contrast does not leave any degree of freedom for the residuals",
                            "in the test of the Space factor.\n\n")
                        skip.Manova <- TRUE
                    }
                } 
                else {
                    dbMEM.S.tmp <- dbmem(sitesX, MEM.autocor="positive")
                    SS <- as.matrix(dbMEM.S.tmp)
                    dbMEM.S.thresh <- give.thresh(dist(sitesX))
                    if(print.res) cat(" Truncation level for space dbMEMs =",
                                      dbMEM.S.thresh, "\n")
                    dbMEM.S <- SS
                    for(j in 2:tt) dbMEM.S <- rbind(dbMEM.S,SS)
                    if(nS==-1) nS <- ncol(SS)
                    else {
                        if(nS > ncol(SS)) {
                            cat("\n Number of requested spatial variables nS =", nS, 
                                "larger than available\n")
                            cat(" The", ncol(SS), 
                                "available spatial dbMEMs will be used\n\n")
                            nS <- ncol(SS)
                        }
                        else dbMEM.S <- dbMEM.S[, 1:nS, drop = FALSE]
                    }
                }
            } else {
                dbMEM.S <- apply(as.matrix(COD.S),2,scale,center=TRUE,scale=TRUE)
                nS <- dim(dbMEM.S)[2]
                if(nS>=s) 
                    stop("The number of spatial coding functions must be lower than S",
                         call.=FALSE) 
                if(nrow(as.matrix(COD.S))!=dim(Y)[1]) 
                    stop("The number of rows in COD.S must be equal to the number",
                         "of observations", call.=FALSE) 
            }
        }
        
        # Generate dbMEM variables for Time (if necessary)
        if(model=="3b" || model=="6b"|| model=="7"||model=="4" || model=="5") {
            if(is.null(COD.T)) {  # Generate temporal dbMEM variables if not given by user
                if(tt==2) {  # Two points in time: construct 2 Helmert contrast variables
                    dbMEM.T=as.matrix(c(rep(-1, s), rep(1,s)))  #DB
                    #
                    if(model=="6b") { 
                        cat("Model 6b requires that 'tt' be larger than 2. When tt=2, full",
                            "coding of the times by a binary variable or Helmert",
                            "contrast does not leave any degree of freedom for the residuals",
                            "in the test of the Time factor.\n\n")
                        skip.Manova <- TRUE
                    }
                    
                } else {  # COD.T is NULL, tt>2
                    if(print.res) cat(" Computing dbMEMs to code for time\n")
                    TT <- as.matrix(dbmem(timesX))
                    dbMEM.T.thresh <- give.thresh(dist(timesX))
                    if(print.res) cat(" Truncation level for time dbMEMs =", 
                                      dbMEM.T.thresh, "\n\n")			      
                    
                    beg <- seq(1,n,by=s)    #  1  7 13 19 25
                    dbMEM.T <- matrix(0,n,ncol(TT))     # Empty matrix 30x2 to house dbMEM.T		
                    for(i in 1:tt) { # Fill tt blocks of rows with identical MEM.T values
                        for(j in 1:s) dbMEM.T[(beg[i]+j-1),] <- TT[i,]
                    }
                    
                    if(nT==-1) { nT <- ncol(TT)   # nT: number of dbMEM variables in TT
                    } else {          # If the arguments specify an imposed value to nT
                        if(nT > ncol(TT)) {
                            cat("\n Number of requested temporal variables nT =", nT, 
                                "larger than available\n")
                            cat(" The", ncol(TT), "available temporal dbMEM will be used\n\n")
                            nT <- ncol(TT)
                        } else {
                            dbMEM.T <- dbMEM.T[, 1:nT, drop = FALSE]
                        }
                    }			      
                }   # End of  if(is.null(COD.T))
                
            } else {   # COD.T is present 
                dbMEM.T=apply(as.matrix(COD.T),2,scale,center=TRUE,scale=TRUE)
                nT <- dim(dbMEM.T)[2]
                if(nT>=tt) stop("The number of temporal coding functions",
                                "must be lower than Ti", call.=FALSE)
                if(nrow(as.matrix(COD.T))!=dim(Y)[1]) 
                    stop("The number of rows in COD.T must be equal to the number",
                         "of observations", call.=FALSE) 
            }
        }
        
        if(!skip.Manova) {   # Begin of "(!skip.Manova)"
            
            if(print.res) {
                if(model!="2" && model!="3b" && model!="6b") { 
                    if(s==2)   #DB To avoid number = -1
                        cat(" Number of space coding functions = 1 \n\n")
                    else
                        cat(" Number of space coding functions =", nS,'\n\n')
                }
                if(model!="2" && model!="3a" && model!="6a") {
                    if(tt==2)  #DB To avoid number = -1
                        cat(" Number of time coding functions = 1 \n\n")
                    else
                        cat(" Number of time coding functions =", nT, "\n\n")}            		
            }
            
            test.STI <- TRUE
            test.S <- TRUE
            test.T <- TRUE
            
            # Generate matrices X (variables for the factor of interest) 
            # and W (covariables) for each test to be done
            if(model=="5") {
                XSTI <- dbMEM.S*dbMEM.T[,1]
                if(dim(dbMEM.T)[2]>1) 
                    for(j in 2:dim(dbMEM.T)[2]) XSTI <- cbind(XSTI,dbMEM.S*dbMEM.T[,j])
                XS <- HM.S
                XT <- HM.T
                if(print.res) {
                    cat(" MODEL V: HELMERT CONTRAST FOR TESTING MAIN FACTORS.\n")
                    cat("          SPACE AND TIME dbMEMs FOR TESTING INTERACTION.\n")
                }					
            } else if(model=="4") {
                XSTI <- dbMEM.S*dbMEM.T[,1]
                if(dim(dbMEM.T)[2]>1) 
                    for(j in 2:dim(dbMEM.T)[2]) XSTI = cbind(XSTI,dbMEM.S*dbMEM.T[,j])
                XS <- dbMEM.S
                XT <- dbMEM.T
                if(print.res) {
                    cat(" MODEL IV: dbMEMs FOR BOTH SPACE AND TIME.\n")
                }		
            } else if(model=="3a") {
                XSTI <- dbMEM.S*HM.T[,1]
                if(tt>2) for(j in 2:(tt-1)) XSTI <- cbind(XSTI,dbMEM.S*HM.T[,j])
                XS <- dbMEM.S
                XT <- HM.T
                if(print.res) {
                    cat(" MODEL IIIa: dbMEMs FOR SPACE AND HELMERT CONTRASTS FOR TIME.\n")
                }	
            } else if(model=="3b") {
                XSTI <- dbMEM.T*HM.S[,1]
                if(dim(HM.S)[2]>1)
                    for(j in 2:(s-1)) XSTI <- cbind(XSTI,dbMEM.T*HM.S[,j])
                XS <- HM.S
                XT <- dbMEM.T
                if(print.res) {
                    cat(" MODEL IIIb: HELMERT CONTRASTS FOR SPACE AND dbMEMs FOR TIME.\n")
                }		
            } else if(model=="7") {
                XSTI <- HM.S*HM.T[,1]
                if(dim(HM.T)[2]>1) for(j in 2:dim(HM.T)[2]) XSTI <- cbind(XSTI,HM.S*HM.T[,j])
                XS <- dbMEM.S
                XT <- dbMEM.T
                if(print.res) {
                    cat(" MODEL VII: dbMEMs FOR BOTH SPACE AND TIME,",
                        "BUT HELMERT CONTRAST FOR INTERACTION.\n")
                }		
                
            } else if(model=="6a") {
                # BEGIN: New code (PL) to assemble a block diagonal matrix for XS
                XS <- matrix(0,n,tt*nS)       # Empty matrix
                beg.x <- seq(1,n,by=s)        # Beginning rows, vertical
                beg.y <- seq(1,tt*nS,by=nS)   # Beginning columns, horizontal
                tmp <- dbMEM.S[1:s,]
                for(j in 1:tt) 
                    XS[(beg.x[j]):(beg.x[j]+s-1),(beg.y[j]):(beg.y[j]+nS-1)] <- tmp
                
                XT <- HM.T
                XSTI = NULL
                if(print.res) {
                    cat(" MODEL VIa: NESTED MODEL.\n")
                    cat("            TESTING FOR THE EXISTENCE",
                        "OF SPATIAL STRUCTURE (COMMON OR SEPARATE)\n")
                }	
                test.STI <- FALSE
                test.T <- FALSE
                
            } else if(model=="6b") {
                # BEGIN: New code (PL)) to assemble a block diagonal matrix for XT
                beg.x <- seq(1,n,by=s)      # Rows to fill at the beginning
                beg.y <- seq(1,s*nT,by=nT)  # Columns to fill at the beginning
                TT = as.matrix(dbMEM.T[beg.x,])   # dbMEM for times, only once
                XT <- matrix(0,n,s*nT)      # Empty matrix to house the diagonal blocks		
                for(j in 1:s) {      # Write groups of s columns in XT
                    for(i in 1:tt) { # Fill tt rows with T1 to T12 in each group of 3 columns
                        XT[(beg.x[i]+j-1), (beg.y[j]):(beg.y[j]+nT-1)] <- TT[i,1:nT]
                    }
                }   ## End of New code
                
                XS <- HM.S						
                XSTI = NULL
                if(print.res) {
                    cat(" MODEL VIb: NESTED MODEL.\n")
                    cat("            TESTING FOR THE EXISTENCE",
                        "OF TEMPORAL STRUCTURE (COMMON OR SEPARATE).\n")
                }	
                test.STI <- FALSE
                test.S <- FALSE
                
            } else if(model=="2") {
                XS <- HM.S
                XT <- HM.T
                XSTI <- NULL
                if(print.res) {
                    cat(" MODEL II: HELMERT CONTRAST FOR SPACE AND TIME.",
                        "NO INTERACTION TERM.\n")
                }	
                test.STI <- FALSE
            } 		
            if(print.res) {
                cat("   Number of space variables =", dim(XS)[2],'\n')
                cat("   Number of time variables =", dim(XT)[2],'\n')
                if(test.STI) {
                    cat("   Number of interaction variables =", dim(XSTI)[2],'\n')
                    cat("   Number of residual degrees of freedom =", 
                        (s*tt-dim(XS)[2]-dim(XT)[2]-dim(XSTI)[2]-1),"\n\n")
                    if((s*tt-dim(XS)[2]-dim(XT)[2]-dim(XSTI)[2]-1) <= 0) 
                        cat("Not enough residual degrees of freedom for testing.\n")
                } else {
                    cat("   Number of residual degrees of freedom =", 
                        (s*tt-dim(XS)[2]-dim(XT)[2]-1),"\n\n")			
                    if((s*tt-dim(XS)[2]-dim(XT)[2]-1) <= 0) 
                        cat("Not enough residual degrees of freedom for testing.\n")
                }
            }		
            
            # Call Manova by RDA, function manovRDa.R
            res <- manovRDa(Y=Y,s=s,tt=tt,S.mat=XS,T.mat=XT,STI.mat=XSTI, Sfixed=Sfixed, 
                            Tfixed=Tfixed, S.test=test.S, T.test=test.T, STI.test=test.STI,  
                            model = model, nperm=nperm, save.X=save.X)
            
            if(test.STI==TRUE) {
                cat(' Interaction test:   R2 =', round(res$testSTI$R2, 4),
                    '  F =',round(res$testSTI$F, 4),'  P(',nperm,'perm) =',res$testSTI$Prob,'\n')
            }
            if(test.S==TRUE) {
                cat(' Space test:         R2 =', round(res$testS$R2, 4),
                    '  F =',round(res$testS$F, 4),'  P(',nperm,'perm) =',res$testS$Prob,'\n')
            }
            if(test.T==TRUE) {
                cat(' Time test:          R2 =', round(res$testT$R2, 4),
                    '  F =',round(res$testT$F, 4),'  P(',nperm,'perm) =',res$testT$Prob,'\n\n')
            }
            
        }  # End of	"if(!skip.Manova)"
        
    })
    aa[3] <- sprintf("%2f",aa[3])
    if(print.res) {
        cat("-------------------------------------------------------\n")
        cat("         Time for computation =",aa[3]," sec",'\n')
        cat("=======================================================\n\n")
    }
    class(res) <- "sti"
    invisible(res)
}




#' @rdname stimodels  
'quicksti' <- function(Y, S, Ti, nperm=999, alpha = 0.05, COD.S=NULL, COD.T=NULL,
                       save.X=FALSE, print.res=TRUE)
{
    if(!is.logical(print.res)) {
        stop("Wrong operator; 'print.res' should be either 'FALSE' or 'TRUE'.", 
             call.=FALSE)
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
        if(dim(Y)[1] != n) stop("The number of rows in species file is not (S x Ti).", 
                                call.=FALSE) 
        
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
        
        # Generates spatial dbMEMs if not provided by user
        if(is.null(COD.S)) {			
            if(print.res) cat(" Computing dbMEMs to code for space\n")
            if(s==2) {  
                dbMEM.S <- as.matrix(rep(c(-1,1),tt))
                nS <- dim(dbMEM.S)[2]
                if(print.res) cat("\n There are only two sites. A vector of Helmert",
                                  "contrasts will be used to represent them\n\n")
            }  
            else {
                SS <- as.matrix(dbmem(sitesX))
                nS <- ncol(SS)
                dbMEM.S.thresh <- give.thresh(dist(sitesX))
                if(print.res) cat(" Truncation level for space dbMEMs =", dbMEM.S.thresh,"\n")
                dbMEM.S <- SS
                for(j in 2:tt) dbMEM.S = rbind(dbMEM.S,SS)
            }
        } else {
            dbMEM.S <- apply(as.matrix(COD.S),2,scale,center=TRUE,scale=TRUE)
            nS <- dim(dbMEM.S)[2]
            if(nS >= s) stop("The number of spatial coding functions must be smaller",
                             "than S.", call.=FALSE) 
            if(nrow(as.matrix(COD.S))!=dim(Y)[1]) stop("The number of rows in COD.S must",
                                                       "be equal to the number of observations.", call.=FALSE) 
        }
        
        # Generate dbMEM variables for time if not provided by user
        if(is.null(COD.T)) {
            nT <- trunc(tt/2)
            if(print.res) cat(" Computing dbMEMs to code for time\n")
            if(tt==2) {  
                dbMEM.T <- as.matrix(c(rep(-1, s), rep(1,s)))  
                if(print.res) cat("\n There are only two time points. A vector of",
                                  "Helmert contrasts will be used to represent them\n\n")
            } else { 
                TT <- as.matrix(dbmem(timesX))
                nT <- ncol(TT)
                dbMEM.T.thresh <- give.thresh(dist(timesX))
                if(print.res) cat(" Truncation level for time dbMEMs =",dbMEM.T.thresh,"\n\n")
                
                
                beg.x <- seq(1,n,by=s)    # New code (PL)
                nT = ncol(TT)             
                dbMEM.T <- matrix(0,n,nT)      # Empty matrix 30x2 to house dbMEM.T		
                for(i in 1:tt) { # Fill tt blocks of rows with identical MEM.T values
                    for(j in 1:s) dbMEM.T[(beg.x[i]+j-1),] <- TT[i,]
                }
                
            }
        } else {
            dbMEM.T <- apply(as.matrix(COD.T),2,scale,center=TRUE,scale=TRUE)
            nT <- dim(dbMEM.T)[2]
            if(nT >= tt) stop("The number of temporal coding functions must be",
                              "lower than Ti.", call.=FALSE)
            if(nrow(as.matrix(COD.T))!=dim(Y)[1]) stop("The number of rows in",
                                                       "COD.T must be equal to the number of observations.", call.=FALSE) 
        }
        
        if(s*tt-s-tt-nT*nS-1<=0) stop("Not enough degrees of freedom for testing",
                                      "interaction.", call.=FALSE)
        
        
        if(print.res) {
            cat(" Number of space coding functions =", nS,'\n')
            cat(" Number of time coding functions =", nT, "\n\n")
        }
        
        # Generates space and time helmert contrasts
        A <- as.factor(rep(1:s,tt))
        B <- rep(1,s)
        for(i in 2:tt) B <- c(B,rep(i,s))
        B <- as.factor(B)
        HM <- model.matrix(~ A + B, contrasts =list(A="contr.helmert", B="contr.helmert"))
        HM.S <- as.matrix(HM[,2:s])
        HM.T <- as.matrix(HM[,(s+1):(s+tt-1)])
        
        res2 <- NULL
        res3 <- NULL
        
        # Test significance for interaction effect 
        #
        # Defines X (variables for the factor of interest) 
        # and W (covariables) for the space-time test
        #
        XSTI <- dbMEM.S*dbMEM.T[,1]
        if(dim(dbMEM.T)[2]>1) for(j in 2:dim(dbMEM.T)[2]) 
            XSTI <- cbind(XSTI,dbMEM.S*dbMEM.T[,j])
        if(print.res) {
            cat("------------------------------------------\n")
            cat(" Testing space-time interaction (model 5)\n")
            cat("------------------------------------------\n\n")
            cat("   Number of space variables =", dim(HM.S)[2],'\n')
            cat("   Number of time variables =", dim(HM.T)[2],'\n')
            cat("   Number of interaction variables =", dim(XSTI)[2],'\n')
            cat("   Number of residual degrees of freedom =", 
                (s*tt-dim(XSTI)[2]-dim(HM.S)[2]-dim(HM.T)[2]-1),"\n\n")
            if((s*tt-dim(XSTI)[2]-dim(HM.S)[2]-dim(HM.T)[2]-1) == 0) 
                stop("Not enough residual degrees of freedom for testing.",
                     "Try with a smaller number of space or time variables.", call.=FALSE)
        }
        
        res <- manovRDa(Y=Y,s=s,tt=tt,S.mat=HM.S,T.mat=HM.T,STI.mat=XSTI, S.test=TRUE, 
                        T.test=TRUE, STI.test=TRUE, model="5", nperm=nperm, save.X=save.X)
        
        cat(' Interaction test:  R2 =', round(res$testSTI$R2, 4),
            '  F =',round(res$testSTI$F,4),'  P(',nperm,'perm) =',res$testSTI$Prob,"\n")
        
        # Begin models 6a and 6b if Prob<=alpha
        if(res$testSTI$Prob<=alpha) {   
            # Model="6a"
            XS <- dbMEM.S
            for(j in 1:(tt-1)) XS = cbind(XS,dbMEM.S*HM.T[,j])
            XT <- HM.T
            if(print.res) {
                cat("----------------------------------------------------\n")
                cat(" Testing for separate spatial structures (model 6a)\n")
                cat("----------------------------------------------------\n\n")
                cat("   Number of space variables =", dim(XS)[2],'\n')
                cat("   Number of time variables =", dim(XT)[2],'\n')
                cat("   Number of residual degrees of freedom =", 
                    (s*tt-dim(XS)[2]-dim(XT)[2]-1),"\n\n")
            }			
            if(s==2) {   # skip.manova
                if(!print.res) cat("\n")
                cat("Model 6a requires that 's' be larger than 2. When s=2, full",
                    "coding of the site locations by a binary variable or Helmert",
                    "contrast does not leave any degree of freedom for the residuals",
                    "in the test of the Space factor.\n\n")
                res2$X.matrix <- NA
                
            } else {
                res2 <- manovRDa(Y=Y,s=s,tt=tt,S.mat=XS,T.mat=XT,STI.mat=NULL, 
                                 S.test=TRUE, T.test=FALSE, STI.test=FALSE, model="6a", 
                                 nperm=nperm, save.X=save.X)
                res$testS <- res2$testS
                cat(' Space test:  R2 =', round(res$testS$R2,4), '  F =',
                    round(res$testS$F,4),'  P(',nperm,'perm) =',res$testS$Prob,"\n")
            }  # End of	"if(s==2)"
            
            # Model="6b"
            XT <- dbMEM.T
            for(j in 1:(s-1)) XT <- cbind(XT,dbMEM.T*HM.S[,j])
            XS <- HM.S
            if(print.res) {
                cat("-----------------------------------------------------\n")
                cat(" Testing for separate temporal structures (model 6b)\n")
                cat("-----------------------------------------------------\n\n")
                cat("   Number of space variables =", dim(XS)[2],'\n')
                cat("   Number of time variables =", dim(XT)[2],'\n')
                cat("   Number of residual degrees of freedom =", 
                    (s*tt-dim(XS)[2]-dim(XT)[2]-1),"\n\n")
            }	
            if(tt==2) {   # skip.manova
                if(!print.res) cat("\n")
                cat("Model 6b requires that 'tt' be larger than 2. When tt=2, full",
                    "coding of the times by a binary variable or Helmert",
                    "contrast does not leave any degree of freedom for the residuals",
                    "in the test of the Time factor.\n\n")
                res3$X.matrix <- NA
            } else {
                res3 <- manovRDa(Y=Y,s=s,tt=tt,S.mat=XS,T.mat=XT,STI.mat=NULL, 
                                 S.test=FALSE, T.test=TRUE, STI.test=FALSE, model="6b", 
                                 nperm=nperm, save.X=save.X)
                res$testT <- res3$testT
                cat(' Time test:   R2 =', round(res$testT$R2,4),'  F =',
                    round(res$testT$F,4),'  P(',nperm,'perm) =',res$testT$Prob,"\n\n")		
            }  # End of	"if(tt==2)"	 	
            
        } else { # End of "if(res$testSTI$Prob<=alpha)",       # End of Models 6a and 6b
            save.X = FALSE
            
            # Tests of Space and Time from the results of model="5"
            if(print.res) {	
                cat(
                    "---------------------------------------------------------------------\n")
                cat(
                    " Testing for common spatial and common temporal structures (model 5)\n")
                cat(
                    "---------------------------------------------------------------------\n\n")
            }
            cat(' Space test:   R2 =', round(res$testS$R2,4),
                '  F =',round(res$testS$F,4),'  P(',nperm,'perm) =',res$testS$Prob,'\n')
            cat(' Time test:    R2 =', round(res$testT$R2,4),
                '  F =',round(res$testT$F,4),'  P(',nperm,'perm) =',res$testT$Prob,'\n\n')
        } 		
    })
    
    aa[3] <- sprintf("%2f",aa[3])
    if(print.res) {
        cat("---------------------------------------------------------\n")
        cat("          Time for computation =",aa[3]," sec",'\n')
        cat("=========================================================\n\n")
    }
    if(save.X) {
        out <- list(res, res2$X.matrix, res3$X.matrix)
        names(out) <- c("tests", "space.blocks", "time.blocks")
        invisible(out)
    }
    else {
        invisible(res)
    } 
}




'manovRDa' <- function(Y, s, tt, S.mat=NULL, T.mat=NULL, STI.mat=NULL, Sfixed=TRUE, 
                       Tfixed=TRUE, S.test=TRUE, T.test=TRUE, STI.test=TRUE, model="5", 
                       nperm=999, save.X)
{
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
    } else {
        projAxB=NULL  
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
    
    # If save.X is TRUE, save the block-diagonal matrix X used in test of model 6a or 6b
    X <- NA
    if(save.X && model=="6a") X <- S.mat
    if(save.X && model=="6b") X <- T.mat
    
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
    
    # Test interaction (unrestricted permutations)
    if(STI.test==TRUE) { 
        nPGE.AxB = 1
        Fref.AxB <- MS.AxB/MS.Res
        
        nPGE.AxB=
            .Call("sti_loop",nperm,Y,s,tt,a,b,cc,SS.Y,Fref.AxB,projAxB,projABAxB)  ### NM
        P.AxB <- nPGE.AxB/(nperm+1)
        
        R2 <- SS.Yfit.AxB/SS.Y
        R2a <- 1-((n-1)/(n-dim(STI.mat)[2]-1))*(1-R2)
        
        testSTI <- 
            list(MS.num=MS.AxB, MS.den=MS.Res, R2=R2, R2.adj=R2a, F=Fref.AxB, Prob=P.AxB)
    } else {
        testSTI <- NULL
    }
    
    
    # Test factor A (space) using restricted permutations within time blocks
    if(S.test==TRUE) { 
        
        #########################
        nPGE.A=1
        if(Tfixed==FALSE && !is.null(STI.mat)) { 
            # Time random factor in crossed design with interaction
            Fref.A <- MS.A/MS.AxB
            MS.den <- MS.AxB
            T_fixed <- 1
        } else if(Tfixed==FALSE && model=="6b") { 
            # Time random factor in nested design
            Fref.A <- MS.A/MS.B
            MS.den <- MS.B
            T_fixed <- 2
        } else {
            Fref.A <- MS.A/MS.Res
            MS.den <-  MS.Res
            T_fixed <- 3
        }
        
        nPGE.A=.Call("s_loop",nperm,Y,s,tt,a,b,cc,SS.Y,Fref.A,projA,projB,projAxB,
                     projABAxB,T_fixed)  ### NM
        
        P.A <- nPGE.A/(nperm+1)
        
        R2 <- SS.Yfit.A/SS.Y
        R2a <- 1-((n-1)/(n-a-1))*(1-R2)
        
        testS <- list(MS.num=MS.A, MS.den=MS.den, R2=R2, R2.adj=R2a, F=Fref.A, Prob=P.A)
    } else {
        testS <- NULL
    }
    
    # Test factor B (time) using restricted permutations within time blocks
    if(T.test==TRUE) { 
        nPGE.B = 1
        if(Sfixed==FALSE && !is.null(STI.mat)) { 
            # Space random factor in crossed design with interaction
            Fref.B <- MS.B/MS.AxB
            MS.den <- MS.AxB
            T_fixed=1   
        } else if(Sfixed==FALSE && model=="6a") { # Space random factor in nested design
            Fref.B <- MS.B/MS.A
            MS.den <- MS.A
            T_fixed=2    
        } else {
            Fref.B <- MS.B/MS.Res
            MS.den <- MS.Res
            T_fixed=3     
        }
        
        
        nPGE.B <- .Call("t_loop",nperm,Y,s,tt,a,b,cc,SS.Y,Fref.B,projA,projB,projAxB,
                        projABAxB,T_fixed)  # NM
        
        P.B <- nPGE.B/(nperm+1)
        
        R2 <- SS.Yfit.B/SS.Y
        R2a <- 1-((n-1)/(n-b-1))*(1-R2)
        
        testT <- list(MS.num=MS.B, MS.den=MS.den, R2=R2, R2.adj=R2a, F=Fref.B, Prob=P.B)
    } else {
        testT <- NULL
    }
    
    return(list(testSTI = testSTI, testS = testS, testT=testT, X.matrix=X))
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
    # restPerm == 2: Restricted permutation of the observations across blocks (within the  
    # blocks formed by each nth observation).
    #
# Vector 'vec' contains the initial order of the objects, e.g. vec=c(1:n).
# At the end of the function, it gives the permuted order of the objects.
#
# Examples:  toto0 <- restrictedPerm(6,4,24,0,c(1:24))
#            toto1 <- restrictedPerm(6,4,24,1,c(1:24))
#
#            Pierre Legendre, January 2006
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
