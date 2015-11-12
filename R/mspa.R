#' Multi-Scale Pattern Analysis
#' 
#' The multi-scale pattern analysis (MSPA, Jombart et al 2009) investigates the 
#' main scales of spatial variation in a multivariate dataset. This 
#' implementation allows one to perform a MSPA using any multivariate analysis 
#' (stored as a \code{\link[ade4]{dudi}} object), and a list of spatial weights
#' (class \code{listw}) or an object of class \code{orthobasisSp}.
#' 
#' The \code{scatter} method is used for plotting the results. Compared to the 
#' original version of the method, this new implementation allows to specify a 
#' number of blocks (\code{nblocks}). In this case, the multiscale decomposition
#' is performed by dividing MEMs into several blocks and summing R2 values. This
#' could facilitate the interpretation of results.
#' 
#' @aliases mspa  print.mspa scatter.mspa
#' @param dudi a duality diagram (i.e. a reduced space ordination) obtained by a
#'   \code{\link[ade4]{dudi}} function (for instance 
#'   \code{\link[ade4]{dudi.pca}}).
#' @param lwORorthobasisSp either a list of weights (class \code{listw}) that 
#'   san be obtained easily using the function \code{\link{chooseCN}} or an 
#'   object of class \code{orthobasisSp}
#' @param nblocks an integer indicating the number of blocks to divide MEMs.
#' @param scannf logical, indicating whether the screeplot should be displayed 
#'   to choose the number or retained factors.
#' @param nf the number of retained factors
#' @param centring a character string indicating if parametric ("param") or 
#'   non-parametric ("sim") centring should be used
#' @param nperm an integer giving the number of permutations used to compute the
#'   theoretical coefficients of determination (999 by default); used if 
#'   centring="sim".
#' @param x a mspa object.
#' @param xax an integer indicating the x axis to be displayed.
#' @param yax an integer indicating the y axis to be displayed.
#' @param posieig a character indicating the position of the screeplot (any of 
#'   the four combination between "top", "bottom", "left" and "right").
#' @param bary a logical indicating whether the barycenter of the variables 
#'   should be displayed.
#' @param plot a logical indicating if the graphics is displayed
#' @param storeData a logical indicating if the data should be stored in the
#'   returned object. If \code{FALSE}, only the names of the data arguments are
#'   stored
#' @param pos an integer indicating the position of the environment where the
#'   data are stored, relative to the environment where the function is called.
#'   Useful only if \code{storeData} is \code{FALSE}
#' @param \dots additional graphical parameters (see ‘adegpar’ and 
#'   ‘trellis.par.get’)
#' @return An object having the classes \code{mspa} and 
#'   \code{\link[ade4]{dudi}}: \code{mspa} objects are \code{\link[ade4]{dudi}}
#'   objects with the following extra slots:\cr - ls: principal components of
#'   the MSPA. These are the coordinates of variables onto principal axes, to be
#'   used for plotting. Correspond to matrix \bold{B} in Appendix A of Jombart
#'   et al (2009). \cr - R2: matrix of R2 between variables and MEMs.
#'   Corresponds to \bold{S} in Jombart et al (2009).\cr - meanPoint:
#'   coordinates of the 'mean variable' onto principal axes. The 'mean variable'
#'   is an hypothetic variable whose scale profile is the average of those of
#'   all variables of the analysed dataset.\cr - varweights: the weights of
#'   variables. Corresponds to \bold{d} in Jombart et al. (2009).\cr
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso \code{\link{chooseCN}} to obtain a list of spatial weights.
#' @references Jombart T, Dray S, and Dufour, A-B. (2009) Finding essential 
#'   scales of spatial variation in ecological data: a multivariate approach. 
#'   \emph{Ecography} \bold{32}: 161-168.
#' @keywords multivariate spatial
#' @examples
#' 
#' 
#' ####################################
#' ### using oribatib mites dataset ###
#' ####################################
#' 
#' if(require("ade4", quietly = TRUE)){
#' ## load data
#' data(oribatid)
#' 
#' ## get the list of spatial weights
#' cn <- chooseCN(oribatid$xy, res = "listw", ask = FALSE, type = 1)
#' 
#' ## Hellinger transformation
#' hellTrans <- function(X){
#'   if (!( is.matrix(X) | is.data.frame(X) )) stop("Object is not a matrix.")  
#'   if (any(is.na(X))) stop("na entries in table.")
#'   
#'   sumRow <- apply(X,1,sum)
#'   Y <- X/sumRow
#'   Y <- sqrt(Y)
#'   
#'   return(Y)
#' }
#' 
#' 
#' ## ENVIRONMENTAL VARIABLES ##
#' ## Hill and Smith analysis for environmental variables
#' ## (for a mixture of quantitative / qualitative variables)
#' hsEnv <- dudi.hillsmith(oribatid$envir,scannf=FALSE)
#' 
#' ## detrending of the analysis (residuals of regression onto xy coordinates)
#' hsEnv.detr <- pcaivortho(hsEnv,oribatid$xy,scannf=FALSE)
#' 
#' ## MSPA of the detrended analysis
#' mspaEnv <- mspa(hsEnv.detr,cn,scannf=FALSE,nf=2)
#' scatter(mspaEnv)
#' 
#' 
#' 
#' ## SPECIES DATA ##
#' ## PCA of species abundances, after Hellinger transformation
#' pcaFau <- dudi.pca(hellTrans(oribatid$fau),scale=FALSE,scannf=FALSE)
#' 
#' ## detrending of this PCA
#' pcaFau.detr <- pcaivortho(pcaFau,oribatid$xy,scannf=FALSE)
#' 
#' # MSPA of the detrended analysis
#' mspaFau <- mspa(pcaFau.detr,cn,scannf=FALSE,nf=2)
#' scatter(mspaFau)
#' 
#' 
#' 
#' ## CANONICAL MSPA ##
#' ## RDA species ~ envir
#' ## (species abundances predicted by environment)
#' ## note: RDA = 'PCAIV' (PCA with Instrumental Variables)
#' rda1 <- pcaiv(dudi=pcaFau.detr, df=oribatid$envir,scannf=FALSE,nf=2)
#' 
#' ## canonical MSPA (species predicted by environment)
#' mspaCan1 <- mspa(dudi=rda1, lw=cn, scannf=FALSE, nf=2)
#' scatter(mspaCan1)
#' 
#' ## same analysis, using a non-parametric centring
#' mspaCan1NP <- mspa(dudi=rda1, lw=cn, scannf=FALSE, nf=2,cent="sim",nper=999)
#' scatter(mspaCan1NP) # basically no change
#' 
#' 
#' 
#' ## PARTIAL CANONICAL MSPA ##
#' ## partial RDA species ~ envir
#' ## (species abundances not predicted by environment)
#' rda2 <- pcaivortho(dudi=pcaFau.detr,df=oribatid$envir,scannf=FALSE,nf=2)
#' 
#' ## partial canonical MSPA
#' mspaCan2 <- mspa(dudi=rda2, lw=cn, scannf=FALSE, nf=2)
#' scatter(mspaCan2) # nothing left
#' }
#' @importFrom ade4 scatter
#' @importFrom ade4 as.dudi scalewt
#' @importFrom lattice xyplot
#' @importFrom stats weighted.mean
#' @importFrom adegraphics sortparamADEgS s.arrow s.label plotEig s1d.barchart
#' @importMethodsFrom adegraphics superpose insert
#' @export mspa

mspa <- function(dudi, lwORorthobasisSp, nblocks, scannf = TRUE, nf = 2, centring = c("param", "sim"), nperm = 999){
    
    ## arguments checks
    if (!inherits(dudi, "dudi"))
        stop("object of class 'dudi' expected")
    if(inherits(lwORorthobasisSp,"orthobasisSp")){
        if (any((dudi$lw - attr(lwORorthobasisSp,"weights"))^2 > 1e-07)) 
            stop("Non equal row weights")
        U <- lwORorthobasisSp
    }
    
    if(inherits(lwORorthobasisSp,"listw"))
        U <- scores.listw(lwORorthobasisSp, wt = dudi$lw, MEM.autocor = "all")
    
    if(ncol(U) != (nrow(U) - 1))
        stop(paste("The orthobasis contains only", ncol(U), "vectors. The decomposition of variance is thus incomplete."))
    
    df <- dudi$tab
    varweights <- dudi$cw/sum(dudi$cw)
    n <- nrow(df)
    p <- ncol(df)
    centring <- match.arg(centring)
    
    ## retrieve var.idx with correct variable names
    findname <- function(vec){
        if(length(vec) == 1) return(vec)
        res <- sub("[.][^.]*$","",vec[1])
        return(res)
    }
    
    ## var.idx
    var.idx <- dudi$assign
    if(!is.null(var.idx)){
        temp <- split(colnames(df), var.idx)
        newlev <- sapply(temp, findname)
        levels(var.idx) <- newlev
    }
    
    # matrix centring and scaling
    X <- scalewt(df, wt = dudi$lw, center=TRUE, scale=TRUE)
    
    ## projection of X onto U
    ## only R-squared of each vector of U are kept
    
    ## model computations
    R <- t(X) %*% diag(dudi$lw) %*% as.matrix(U)
    R2 <- R*R
    
    ## handle centring of R2
    if(centring=="param"){
        newdf <- R2-(1/(n-1))
    } else{ # i.e. if centring is non-parametric
        tempX <- t(X)
        fPerm <- function(X){
            permX <- X[, sample(1:n)] # X has to be transposed, that is, variables in rows, obs in columns
            res <- permX %*% diag(dudi$lw) %*% as.matrix(U) 
            res <- res * res
            return(res)
        } # end fPerm
        listR2sim <- lapply(1:nperm, function(i) fPerm(tempX))
        meanR2sim <- listR2sim[[1]]
        for(i in 2:nperm){
            meanR2sim <- meanR2sim + listR2sim[[i]]
        } # end for
        
        meanR2sim <- meanR2sim / nperm
        
        newdf <- R2 - meanR2sim
    } # end non-parametric centring
    
    ## keep only positive deviation in R2
    newdf[newdf < 0] <- 0
    
    ## smoothed analysis (using blocks of MEMs)
    if(!(missing(nblocks))){
        fac <- cut(1:ncol(U), nblocks)
        i.start <- tapply(1:ncol(U), fac, min)
        i.stop <- tapply(1:ncol(U), fac, max)
        levels(fac) <- paste("[", i.start, "-", i.stop, "]", sep="")
        newdf <- t(apply(newdf, 1, function(i) tapply(i, fac, sum)))
        R2 <- t(apply(R2, 1, function(i) tapply(i, fac, sum)))
    }
    newdf <- as.data.frame(newdf)
    
    ## we proceed to the analysis of this matrix
    res <- as.dudi(newdf, scannf = scannf, nf = nf, row.w = varweights,
                   col.w = rep(1, ncol(newdf)), call = match.call(), type = "mspa")
    
    res$ls <- as.data.frame(as.matrix(R2) %*% as.matrix(res$c1))
    colnames(res$ls) <- colnames(res$li)
    row.names(res$ls) <- row.names(res$li)
    
    xmoy <- apply(R2, 2, function(c) weighted.mean(c,varweights))
    bary <- as.vector(t(xmoy) %*% as.matrix(res$c1))
    names(bary) <- colnames(res$c1)
    
    res$R2 <- R2
    res$meanPoint <- bary
    res$varweights <- varweights
    names(res$varweights) <- colnames(X)
    if(!is.null(var.idx)) res$assign <- var.idx
    if(centring=="sim") {
        res$centring <- meanR2sim
        if(!(missing(nblocks)))
            res$centring <- tapply(meanR2sim, fac, sum)
    }
    
    return(res)
} 

#' @rdname mspa
#' @export
scatter.mspa <- function(x, xax = 1, yax = 2, posieig = "topleft", bary=TRUE, 
                         plot = TRUE, storeData = TRUE, pos = -1, ...){
    
    if(!inherits(x,"mspa"))
        stop("Object of class 'mspa' expected")
    
    if((xax == yax) || (x$nf == 1))
        stop("One axis only : not yet implemented")
    if(length(xax) > 1 | length(yax) > 1)
        stop("Not implemented for multiple xax/yax")
    
    if(xax > x$nf)
        stop("Non convenient xax")
    if(yax > x$nf)
        stop("Non convenient yax")  
    
    position <- match.arg(posieig[1], choices = c("bottomleft", "bottomright", "topleft", "topright", "none"), several.ok = FALSE)
    
    ## sort parameters for each graph
    graphsnames <- c("mem", "var", "bary", "eig")
    sortparameters <- sortparamADEgS(..., graphsnames = graphsnames)
    
    ## parameters management
    params <- list()
    params$mem <- list(plabels = list(cex = 1))
    params$var <- list(plabels = list(cex = 1.5))
    params$eig <- list(pbackground = list(box = TRUE), psub = list(text = "Eigenvalues"))
    if(bary)
        params$bary = list(ppoints = list(pch = 20, cex = 2, col = "black"))
    sortparameters <- modifyList(params, sortparameters, keep.null = TRUE)
    
    g1 <- do.call("s.arrow", c(list(dfxy = substitute(x$c1), xax = xax, yax = yax, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$mem))
    g2 <- do.call("s.label", c(list(dfxy = substitute(x$ls), xax = xax, yax = yax, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$var))
    object <- do.call("superpose", list(g1@Call, g2@Call))
    
    if(bary){
        g3 <- xyplot(x$meanPoint[yax] ~ x$meanPoint[xax], col = sortparameters$bary$ppoints$col, pch = sortparameters$bary$ppoints$pch, cex = sortparameters$bary$ppoints$cex, xlab = "", ylab = "", aspect = g1@adeg.par$paxes$aspectratio)
        object <- do.call("superpose", list(object@Call, g3))
    }
    
    if(position != "none") {
        g4 <- do.call("plotEig", c(list(eigvalue = substitute(x$eig), nf = 1:x$nf, xax = xax, yax = yax, plot = FALSE), sortparameters$eig))
        object <- do.call("insert", list(g4@Call, object@Call, posi = position, ratio = 0.25, plot = FALSE))        
    }
    
    names(object) <- graphsnames[c(1, 2, 3 * isTRUE(bary), 4 * (position != "none"))]
    object@Call <- match.call()
    if(plot) 
        print(object)
    invisible(object)
} 

#' @rdname mspa
#' @export
print.mspa <- function(x, ...){
    cat("Multi-Scale Pattern Analysis\n")
    cat("call: ")
    print(x$call)
    cat("class: ")
    cat(class(x), "\n")
    cat("\n$rank (rank)     :", x$rank)
    cat("\n$nf (axis saved) :", x$nf)
    cat("\n\neigen values: ")
    l0 <- length(x$eig)
    cat(signif(x$eig, 4)[1:(min(5, l0))])
    if (l0 > 5)
        cat(" ...\n\n")
    else cat("\n\n")
    sumry <- array("", c(3, 4), list(1:3, c("vector", "length",
                                            "mode", "content")))
    sumry[1, ] <- c("$eig", length(x$eig), mode(x$eig), "Eigenvalues")
    sumry[2, ] <- c("$cw", length(x$cw), mode(x$cw), "MEM weights")
    sumry[3, ] <- c("$lw", length(x$lw), mode(x$lw), "variables weights")
    
    print(sumry, quote = FALSE)
    cat("\n")
    
    sumryA <- array("", c(4, 4), list(1:4, c("data.frame", "nrow", "ncol", "content")))
    sumryA[1, ] <- c("$R2", nrow(x$R2), ncol(x$R2), "Matrix of R2 ('S')")
    sumryA[2, ] <- c("$tab", nrow(x$tab), ncol(x$tab), "Matrix of centred R2 ('Z')")
    sumryA[3, ] <- c("$c1", nrow(x$c1), ncol(x$c1), "Principal axes, i.e. MEM loadings('A')")
    sumryA[4, ] <- c("$ls", nrow(x$ls), ncol(x$ls), "Projection of variables (R2) on principal axes ('B')")
    
    cat("\nMain components:\n")
    print(sumryA, quote = FALSE)
    
    sumryB <- array("", c(3, 4), list(1:3, c("data.frame", "nrow", "ncol", "content")))
    sumryB[1, ] <- c("$li", nrow(x$li), ncol(x$li), "Scores for variables")
    sumryB[2, ] <- c("$l1", nrow(x$l1), ncol(x$l1), "Principal components")
    sumryB[3, ] <- c("$co", nrow(x$co), ncol(x$co), "Scores for MEM")
    
    cat("\nGeneric 'dudi' components:\n")
    print(sumryB, quote = FALSE)
    
    cat("Other elements: ")
    if (length(names(x)) > 14)
        cat(names(x)[15:(length(x))], "\n")
    else cat("NULL\n")
} 
