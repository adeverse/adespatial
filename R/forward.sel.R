#' Forward selection with multivariate Y using permutation under reducel model
#' 
#' Performs a forward selection by permutation of residuals under reduced model.
#' Y can be multivariate.
#' 
#' The forward selection will stop when either K, R2tresh, adjR2tresh, alpha and
#' R2more has its parameter reached.
#' 
#' @aliases forward.sel
#' @param Y Response data matrix with n rows and m columns containing 
#'   quantitative variables
#' @param X Explanatory data matrix with n rows and p columns containing 
#'   quantitative variables
#' @param K Maximum number of variables to be selected. The default is one minus
#'   the number of rows
#' @param R2thresh Stop the forward selection procedure if the R-square of the 
#'   model exceeds the stated value. This parameter can vary from 0.001 to 1
#' @param adjR2thresh Stop the forward selection procedure if the adjusted 
#'   R-square of the model exceeds the stated value. This parameter can take any
#'   value (positive or negative) smaller than 1
#' @param R2more Stop the forward selection procedure if the difference in model
#'   R-square with the previous step is lower than R2more. The default setting 
#'   is 0.001
#' @param alpha Significance level. Stop the forward selection procedure if the 
#'   p-value of a variable is higher than alpha. The default is 0.05 is TRUE
#' @param nperm The number of permutation to be used.The default setting is 999 
#'   permutation.
#' @param Xscale Standardize the variables in table X to variance 1. The default
#'   setting is TRUE
#' @param Ycenter Center the variables in table Y. The default setting is TRUE
#' @param Yscale Standardize the variables in table Y to variance 1. The default
#'   setting is FALSE.
#' @param verbose If 'TRUE' more diagnostics are printed. The default setting is
#'   TRUE
#' @return A dataframe with: \item{ variables }{ The names of the variables } 
#'   \item{ order }{ The order of the selection of the variables } \item{ R2 }{ 
#'   The R2 of the variable selected } \item{ R2Cum }{ The cumulative R2 of the 
#'   variables selected } \item{ AdjR2Cum }{ The cumulative adjusted R2 of the 
#'   variables selected } \item{ F }{ The F statistic } \item{ pval }{ The 
#'   P-value statistic }
#' @note Not yet implemented for CCA (weighted regression) and with covariables.
#' @author Stephane Dray \email{stephane.dray@@univ-lyon1.fr}
#' @references Canoco manual p.49
#'   
#' @keywords multivariate
#' @examples
#' 
#' x <- matrix(rnorm(30),10,3)
#' y <- matrix(rnorm(50),10,5)
#'     
#' forward.sel(y,x,nperm=99, alpha = 0.5)
#'  
#' @useDynLib adespatial, .registration = TRUE 
#' @export forward.sel
"forward.sel" <-
    function(Y, X , K = nrow(X)-1, R2thresh = .99, adjR2thresh = .99, nperm = 999, R2more = 0.001,
        alpha = 0.05, Xscale = TRUE, Ycenter = TRUE, Yscale = FALSE, verbose = TRUE){
        X <- as.data.frame(X)
        Y <- as.data.frame(Y)
        if (any(is.na(X))|any(is.na(X))) stop("na entries in table")
        if(nrow(X)!=nrow(Y)) stop("different number of rows")
        if (any(apply(X,2,is.factor))|any(apply(Y,2,is.factor))) stop("not yet implemented for factors")
        X <- apply(X,2,scale,scale=Xscale)
        Y <- apply(Y,2,scale,scale=Yscale,center=Ycenter)
        nbcovar <- 0
        ##W <- NULL
        ##if(!is.null(W)){
        ##  W <- as.data.frame(W)
        ##  W <- apply(W,2,scale,scale=Wscale)
        ##  if(nrow(W)!=nrow(Y)) stop("different number of rows")
        ##  Yori <- Y
        ##  Xori <- X
        ##  X <- as.data.frame(residuals(lm(as.matrix(X)~as.matrix(W))))
        ##  Y <- as.data.frame(residuals(lm(as.matrix(Y)~as.matrix(W))))
        ##  nbcovar=ncol(W)
        ##}
        pval <- rep(1,ncol(X))
        ordre <- rep(0,ncol(X))
        R2 <- rep(0,ncol(X))
        adjR2 <- rep(0,ncol(X))
        Fvalue <- rep(0,ncol(X))
        res <- list()
        res <- .C("forwardsel", as.double(t(X)), as.double(t(Y)), as.integer(nrow(X)),as.integer(ncol(X)),as.integer(ncol(Y)),pval=as.double(pval), ord=as.integer(ordre),Fval=as.double(Fvalue),as.integer(nperm),R2=as.double(R2),adjR2=as.double(adjR2),as.integer(K),as.double(R2thresh),as.double(adjR2thresh),as.double(R2more),as.integer(nbcovar),as.double(alpha), as.integer(verbose), PACKAGE="adespatial")[c("ord","Fval","pval","R2","adjR2")]
        lambdA <- c(res$R2[1],diff(res$R2))
        resmat <- data.frame(res$ord,lambdA,res$R2,res$adjR2,res$Fval,res$pval)
        if(sum(res$ord>0)==0) stop("No variables selected. Please change your parameters.")
        resmat <- resmat[res$ord>0,]
        resmat <- cbind(I(colnames(X)[resmat[,1]]),resmat)
        ##if(!is.null(W)){
        ##  Yori <- as.matrix(Yori)
        ##  Y <- as.matrix(Y)
        ##  trori <- sum(diag(crossprod(Yori)))
        ##  trdt <- sum(diag(crossprod(Y)))
        ##  resmat[,3] <- resmat[,3]*trdt/trori
        ##  resmat[,4] <- resmat[,4]*trdt/trori
        ##}
        names(resmat) <- c("variables","order","R2","R2Cum","AdjR2Cum","F","pvalue")
        return(resmat)
    }
