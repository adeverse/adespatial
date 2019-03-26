#' Multivariate spatial analysis
#' 
#' This function provides a multivariate extension of the univariate method of
#' spatial autocorrelation analysis. It provides a spatial ordination by maximizing
#' the product of variance by spatial autocorrelation.
#' 
#' This analysis generalizes the Wartenberg's multivariate spatial correlation
#' analysis to various duality diagrams created by the functions
#' (\code{dudi.pca}, \code{dudi.coa}, \code{dudi.acm}, \code{dudi.mix}...) If
#' \emph{dudi} is a duality diagram created by the function \code{dudi.pca} and
#' \emph{listw} gives spatial weights created by a row normalized coding
#' scheme, the analysis is equivalent to Wartenberg's analysis.
#' 
#' We note X the data frame with the variables, Q the column weights matrix and
#' D the row weights matrix associated to the duality diagram \emph{dudi}. We
#' note L the neighbouring weights matrix associated to \emph{listw}. Then, the
#' \code{'multispati'} analysis gives principal axes v that maximize the
#' product of spatial autocorrelation and inertia of row scores :
#' \deqn{I(XQv)*\|XQv\|^2 = v^{t}Q^{t}X^{t}DLXQv}{I(XQv)*\|\|XQv\|\|^2 =
#' t(v)t(Q)t(X)DLXQv}
#' 
#' @aliases multispati plot.multispati summary.multispati print.multispati
#' @param dudi an object of class \code{dudi} obtained by the simple analysis of a data table
#' @param listw an object of class \code{listw} created for example by \code{\link[spdep]{nb2listw}}
#' @param scannf a logical value indicating whether the eigenvalues barplot
#' should be displayed
#' @param nfposi an integer indicating the number of axes with positive autocorrelation
#' @param nfnega an integer indicating the number of axes with negative autocorrelation
#' @param x,object an object of class \code{multispati}
#' @param xax,yax the numbers of the x-axis and the y-axis
#' @param plot a logical indicating if the graphics is displayed
#' @param storeData a logical indicating if the data should be stored in the
#'   returned object. If \code{FALSE}, only the names of the data arguments are
#'   stored
#' @param pos an integer indicating the position of the environment where the
#'   data are stored, relative to the environment where the function is called.
#'   Useful only if \code{storeData} is \code{FALSE}

#' @param \dots further arguments passed to or from other methods
#' @return Returns an object of class \code{multispati}, which contains the
#' following elements : 
#' \item{eig}{a numeric vector containing the eigenvalues}
#' \item{nfposi}{integer, number of kept axes associated to positive
#' eigenvalues} 
#' \item{nfnega}{integer, number of kept axes associated to
#' negative eigenvalues} 
#' \item{c1}{principle axes (v), data frame with p rows
#' and (nfposi + nfnega) columns} 
#' \item{li}{principal components (XQv), data
#' frame with n rows and (nfposi + nfnega) columns} 
#' \item{ls}{lag vector onto
#' the principal axes (LXQv), data frame with n rows and (nfposi + nfnega)
#' columns} 
#' \item{as}{principal axes of the dudi analysis (u) onto principal
#' axes of multispati (t(u)Qv), data frame with dudi\$nf rows and (nfposi +
#' nfnega) columns}
#' @author Stéphane Dray \email{stephane.dray@@univ-lyon1.fr} with contributions by 
#' Daniel Chessel, Sebastien Ollier and Thibaut Jombart
#' 
#' @seealso \code{\link[ade4]{dudi}},\code{\link[spdep]{mat2listw}}
#' @references Dray, S., Said, S. and Debias, F. (2008) Spatial ordination of
#' vegetation data using a generalization of Wartenberg's multivariate spatial
#' correlation. \emph{Journal of vegetation science}, \bold{19}, 45--56.
#' 
#' Grunsky, E. C. and Agterberg, F. P. (1988) Spatial and multivariate analysis
#' of geochemical data from metavolcanic rocks in the Ben Nevis area, Ontario.
#' \emph{Mathematical Geology}, \bold{20}, 825--861.
#' 
#' Switzer, P. and Green, A.A. (1984) Min/max autocorrelation factors for
#' multivariate spatial imagery. Tech. rep. 6, Stanford University.
#' 
#' Thioulouse, J., Chessel, D. and Champely, S. (1995) Multivariate analysis of
#' spatial patterns: a unified approach to local and global structures.
#' \emph{Environmental and Ecological Statistics}, \bold{2}, 1--14.
#' 
#' Wartenberg, D. E. (1985) Multivariate spatial correlation: a method for
#' exploratory geographical analysis. \emph{Geographical Analysis}, \bold{17},
#' 263--283.
#' 
#' @keywords multivariate spatial
#' @examples
#' 
#' 
#' if (require(spdep, quiet = TRUE) & require(ade4, quiet = TRUE)) {
#'     data(mafragh)
#'     maf.xy <- mafragh$xy
#'     maf.flo <- mafragh$flo
#'     maf.listw <- nb2listw(mafragh$nb)
#'     if(adegraphicsLoaded()) {
#'       g1 <- s.label(maf.xy, nb = mafragh$nb, plab.cex = 0.75)
#'     } else {
#'       s.label(maf.xy, neig = mafragh$neig, clab = 0.75)
#'     }
#'     maf.coa <- dudi.coa(maf.flo,scannf = FALSE)
#'     maf.coa.ms <- multispati(maf.coa, maf.listw, scannf = FALSE, nfposi = 2, nfnega = 2)
#'     maf.coa.ms
#'     
#'     ### detail eigenvalues components
#'     fgraph <- function(obj){
#'       # use multispati summary
#'       sum.obj <- summary(obj)
#'       # compute Imin and Imax
#'       Ibounds <- moran.bounds(eval(as.list(obj$call)$listw))
#'       Imin <- Ibounds[1]
#'       Imax <- Ibounds[2]
#'       I0 <- -1/(nrow(obj$li)-1)
#'       # create labels
#'       labels <- lapply(1:length(obj$eig),function(i) bquote(lambda[.(i)]))
#'       # draw the plot
#'       xmax <- eval(as.list(obj$call)$dudi)$eig[1]*1.1
#'       par(las=1)
#'       var <- sum.obj[,2]
#'       moran <- sum.obj[,3]
#'       plot(x=var,y=moran,type='n',xlab='Inertia',ylab="Spatial autocorrelation (I)",
#'            xlim=c(0,xmax),ylim=c(Imin*1.1,Imax*1.1),yaxt='n')
#'       text(x=var,y=moran,do.call(expression,labels))
#'       ytick <- c(I0,round(seq(Imin,Imax,le=5),1))
#'       ytlab <- as.character(round(seq(Imin,Imax,le=5),1))
#'       ytlab <- c(as.character(round(I0,1)),as.character(round(Imin,1)),
#'            ytlab[2:4],as.character(round(Imax,1)))
#'       axis(side=2,at=ytick,labels=ytlab)
#'       rect(0,Imin,xmax,Imax,lty=2)
#'       segments(0,I0,xmax,I0,lty=2)
#'       abline(v=0)
#'       title("Spatial and inertia components of the eigenvalues")
#'     }
#'     fgraph(maf.coa.ms)
#'     ## end eigenvalues details
#' 
#' 
#'     if(adegraphicsLoaded()) {
#'       g2 <- s1d.barchart(maf.coa$eig, p1d.hori = FALSE, plot = FALSE)
#'       g3 <- s1d.barchart(maf.coa.ms$eig, p1d.hori = FALSE, plot = FALSE) 
#'       g4 <- s.corcircle(maf.coa.ms$as, plot = FALSE)
#'       G1 <- ADEgS(list(g2, g3, g4), layout = c(1, 3))
#'     } else {
#'       par(mfrow = c(1, 3))
#'       barplot(maf.coa$eig)
#'       barplot(maf.coa.ms$eig) 
#'       s.corcircle(maf.coa.ms$as)
#'       par(mfrow = c(1, 1))
#'     }
#'  
#'  
#'     if(adegraphicsLoaded()) {
#'       g5 <- s.value(maf.xy, -maf.coa$li[, 1], plot = FALSE)
#'       g6 <- s.value(maf.xy, -maf.coa$li[, 2], plot = FALSE)
#'       g7 <- s.value(maf.xy, maf.coa.ms$li[, 1], plot = FALSE)
#'       g8 <- s.value(maf.xy, maf.coa.ms$li[, 2], plot = FALSE)
#'       G2 <- ADEgS(list(g5, g6, g7, g8), layout = c(2, 2))
#'     } else {
#'       par(mfrow = c(2, 2))
#'       s.value(maf.xy, -maf.coa$li[, 1])
#'       s.value(maf.xy, -maf.coa$li[, 2])
#'       s.value(maf.xy, maf.coa.ms$li[, 1])
#'       s.value(maf.xy, maf.coa.ms$li[, 2])
#'       par(mfrow = c(1, 1))
#'     }
#' 
#' 
#'     w1 <- -maf.coa$li[, 1:2]
#'     w1m <- apply(w1, 2, lag.listw, x = maf.listw)
#'     w1.ms <- maf.coa.ms$li[, 1:2]
#'     w1.msm <- apply(w1.ms, 2, lag.listw, x = maf.listw)
#'     if(adegraphicsLoaded()) {
#'       g9 <- s.match(w1, w1m, plab.cex = 0.75, plot = FALSE)
#'       g10 <- s.match(w1.ms, w1.msm, plab.cex = 0.75, plot = FALSE)
#'       G3 <- cbindADEg(g9, g10, plot = TRUE)
#'     } else {
#'       par(mfrow = c(1,2))
#'       s.match(w1, w1m, clab = 0.75)
#'       s.match(w1.ms, w1.msm, clab = 0.75)
#'       par(mfrow = c(1, 1))
#'     }
#' 
#'     maf.pca <- dudi.pca(mafragh$env, scannf = FALSE)
#'     multispati.randtest(maf.pca, maf.listw)
#'     maf.pca.ms <- multispati(maf.pca, maf.listw, scannf=FALSE)
#'     plot(maf.pca.ms)
#' }
#' 
#' 
#' @importFrom spdep lag.listw
#' @importFrom adegraphics sortparamADEgS s.match s.arrow plotEig layout2position s.corcircle
#' @importFrom graphics barplot
#' @importFrom methods new
#' @export

"multispati" <- function(dudi, listw, scannf = TRUE, nfposi = 2, nfnega = 0) {
    if(!inherits(dudi,"dudi")) 
        stop ("object of class 'dudi' expected")
    if(!inherits(listw,"listw")) 
        stop ("object of class 'listw' expected") 
    if(listw$style!="W") 
        stop ("object of class 'listw' with style 'W' expected")
    NEARZERO <- 1e-14
    
    dudi$cw <- dudi$cw
    fun <- function (x) lag.listw(listw, x, TRUE)
    tablag <- apply(dudi$tab,2,fun)
    covar <- t(tablag) %*% as.matrix((dudi$tab*dudi$lw))
    covar <- (covar+t(covar))/2
    covar <- covar * sqrt(dudi$cw)
    covar <- t(t(covar) * sqrt(dudi$cw))
    covar <- eigen(covar, symmetric = TRUE)
    res <- list()
    res$eig <- covar$values[abs(covar$values)>NEARZERO]
    ndim <- length(res$eig)
    covar$vectors <- covar$vectors[, abs(covar$values)>NEARZERO]
    
    if (scannf) {
        barplot(res$eig)
        cat("Select the number of axes with positive spatial autocorrelation: ")
        nfposi <- as.integer(readLines(n = 1))
        
        cat("Select the number of axes with negative spatial autocorrelation: ")
        nfnega <- as.integer(readLines(n = 1))
    }
    
    if (nfposi <= 0) nfposi <- 1
    if (nfnega <= 0) nfnega <- 0       
    
    if(nfposi > sum(res$eig > 0)){
          nfposi <- sum(res$eig > 0)
          warning(paste("There are only",sum(res$eig>0),"positive factors."))
        }
    if(nfnega > sum(res$eig < 0)){
          nfnega <- sum(res$eig < 0)
          warning(paste("There are only",sum(res$eig< 0),"negative factors."))
        }
    
    res$nfposi <- nfposi
    res$nfnega <- nfnega
    agarder <- c(1:nfposi, if (nfnega>0) (ndim-nfnega+1):ndim else NULL)
    dudi$cw[which(dudi$cw == 0)] <- 1
    auxi <- data.frame(covar$vectors[, agarder] /sqrt(dudi$cw))
    names(auxi) <- paste("CS", agarder, sep = "")
    row.names(auxi) <- names(dudi$tab)
    res$c1 <- auxi                     
    auxi <- as.matrix(auxi)*dudi$cw
    auxi1 <- as.matrix(dudi$tab)%*%auxi
    auxi1 <- data.frame(auxi1)
    names(auxi1) <- names(res$c1)
    row.names(auxi1) <- row.names(dudi$tab)
    res$li <- auxi1
    auxi1 <- as.matrix(tablag)%*%auxi
    auxi1 <- data.frame(auxi1)
    names(auxi1) <- names(res$c1)
    row.names(auxi1) <-  row.names(dudi$tab)    
    res$ls <- auxi1
    auxi <- as.matrix(res$c1) * unlist(dudi$cw)
    auxi <- data.frame(t(as.matrix(dudi$c1)) %*% auxi)
    row.names(auxi) <- names(dudi$li)
    names(auxi) <- names(res$li)
    res$as <- auxi
    res$call <- match.call()
    class(res) <- "multispati"
    return(res)
}

#' @rdname multispati
#' @export
"summary.multispati" <- function (object, ...) {
  
  norm.w <- function(X, w) {
    f2 <- function(v) sum(v * v * w)/sum(w)
    norm <- apply(X, 2, f2)
    return(norm)
  }
  
  if (!inherits(object, "multispati")) stop("to be used with 'multispati' object")
  
  cat("\nMultivariate Spatial Analysis\n")
  cat("Call: ")
  print(object$call)
  
  appel <- as.list(object$call)
  dudi <- eval.parent(appel$dudi)
  listw <- eval.parent(appel$listw)
  
  ## les scores de l'analyse de base
  nf <- dudi$nf
  eig <- dudi$eig[1:nf]
  cum <- cumsum (dudi$eig) [1:nf]
  ratio <- cum/sum(dudi$eig)
  w <- apply(dudi$l1, 2, spdep::lag.listw, x = listw, zero.policy = TRUE)
  moran <- apply(w*as.matrix(dudi$l1)*dudi$lw,2,sum)
  res <- data.frame(var=eig,cum=cum,ratio=ratio, moran=moran)
  cat("\nScores from the initial duality diagram:\n")
  print(res)
  
  ## les scores de l'analyse spatiale
  ## on recalcule l'objet en gardant tous les axes
  eig <- object$eig
  nfposi <- object$nfposi
  nfnega <- object$nfnega
  nfposimax <- sum(eig > 0)
  nfnegamax <- sum(eig < 0)
  
  ms <- multispati(dudi=dudi, listw=listw, scannf=FALSE,
                   nfposi=nfposimax, nfnega=nfnegamax)
  
  ndim <- dudi$rank
  nf <- nfposi + nfnega
  agarder <- c(1:nfposi,if (nfnega>0) (ndim-nfnega+1):ndim else NULL)
  varspa <- norm.w(ms$li,dudi$lw)
  moran <- apply(as.matrix(ms$li)*as.matrix(ms$ls)*dudi$lw,2,sum)
  res <- data.frame(eig=eig,var=varspa,moran=moran/varspa)
  
  cat("\nMultispati eigenvalues decomposition:\n")
  print(res[agarder,])
  return(invisible(res))
}


#' @rdname multispati
#' @export
print.multispati <- function(x, ...)
{
    cat("Multispati object \n")
    cat("class: ")
    cat(class(x))
    cat("\n$call: ")
    print(x$call)
    cat("\n$nfposi:", x$nfposi, "axis-components saved")
    cat("\n$nfnega:", x$nfnega, "axis-components saved")
    #cat("\n$rank: ")
    #cat(x$rank)
    cat("\nPositive eigenvalues: ")
    l0 <- sum(x$eig >= 0)
    cat(signif(x$eig, 4)[1:(min(5, l0))])
    if (l0 > 5) 
        cat(" ...\n")
    else cat("\n")  
    cat("Negative eigenvalues: ")
    l0 <- sum(x$eig <= 0)
    cat(sort(signif(x$eig, 4))[1:(min(5, l0))])
    if (l0 > 5) 
        cat(" ...\n")
    else cat("\n")
    cat('\n')
    sumry <- array("", c(1, 4), list(1, c("vector", "length", 
        "mode", "content")))
    sumry[1, ] <- c('$eig', length(x$eig), mode(x$eig), 'eigen values')
    
    print(sumry, quote = FALSE)
    cat("\n")
    sumry <- array("", c(4, 4), list(1:4, c("data.frame", "nrow", "ncol", "content")))
    sumry[1, ] <- c("$c1", nrow(x$c1), ncol(x$c1), "column normed scores")
    sumry[2, ] <- c("$li", nrow(x$li), ncol(x$li), "row coordinates")
    sumry[3, ] <- c("$ls", nrow(x$ls), ncol(x$ls), 'lag vector coordinates')
    sumry[4, ] <- c("$as", nrow(x$as), ncol(x$as), 'inertia axes onto multispati axes')
    
    
    print(sumry, quote = FALSE)
    cat("other elements: ")
    if (length(names(x)) > 8) 
        cat(names(x)[9:(length(names(x)))], "\n")
    else cat("NULL\n")
}

#' @rdname multispati
#' @export
"plot.multispati" <- function(x, xax = 1, yax = 2, pos = -1, storeData = TRUE, plot = TRUE, ...) {
    if(!inherits(x, "multispati")) 
        stop("Object of class 'multispati' expected")
    if((xax == yax) || ((x$nfposi + x$nfnega) == 1))
        stop("One axis only : not yet implemented")
    if(length(xax) > 1 | length(yax) > 1)
        stop("Not implemented for multiple xax/yax")
    
    if(xax > (x$nfposi + x$nfnega)) 
        stop("Non convenient xax")
    if(yax > (x$nfposi + x$nfnega)) 
        stop("Non convenient yax")
    
    ## sort parameters for each graph
    graphsnames <- c("row", "eig", "loadings", "Xax")
    sortparameters <- sortparamADEgS(..., graphsnames = graphsnames)
    
    ## default values for parameters
    params <- list()
    params[[1]] <- list(psub = list(text = "Scores and lag scores"))
    params[[2]] <- list(psub = list(text = "Eigenvalues"), paxes = list(draw = TRUE, x = list(draw = FALSE), y = list(draw = TRUE)))
    params[[3]] <- list(psub = list(text = "Loadings"))
    params[[4]] <- list(psub = list(text = "Unconstrained axes"), pbackground = list(box = FALSE), plabels = list(cex = 1.25))
    names(params) <- graphsnames
    sortparameters <- modifyList(params, sortparameters, keep.null = TRUE)
    
    ## creation of each individual ADEg
    g1 <- do.call("s.match", c(list(dfxy1 = substitute(x$li), dfxy2 = substitute(x$ls), xax = xax, yax = yax, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters[[1]]))
    g2 <- do.call("plotEig", c(list(eigvalue = substitute(x$eig), nf = c(1:x$nfposi, length(x$eig):(length(x$eig) - x$nfnega + 1)), xax = xax, yax = yax, plot = FALSE), sortparameters[[2]]))
    g3 <- do.call("s.arrow", c(list(dfxy = substitute(x$c1), xax = xax, yax = yax, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters[[3]]))
    g4 <- do.call("s.corcircle",c(list(dfxy = substitute(x$as), xax = xax, yax = yax, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters[[4]]))
    
    ## ADEgS creation
    lay <-  matrix(c(rep(0, 4), 2, 2, rep(1, 4), 2, 2, rep(1, 4), 3, 3, rep(1, 4), 3, 3, rep(1, 4), 4, 4, rep(0, 4), 4, 4), 6, 6) 
    object <- new(Class = "ADEgS", ADEglist = list(g1, g2, g3, g4), positions = layout2position(lay), add = matrix(0, ncol = 4, nrow = 4), Call = match.call())
    names(object) <- graphsnames
    if(plot)
        print(object)
    invisible(object)   
}

