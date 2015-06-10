WRperiodogram <- function(x, T1 = 2, T2, nperm = 499, nopermute, mult = c("sidak", "bonferroni"), print.time = FALSE) 
{
  if(!is.numeric(x)) stop("x contains non-numeric values; it should only contain real numbers")
  mult <- match.arg(mult)
  n <- length(x)
  T.sup <- n %/% 2
  if(missing(T2)) T2 <- T.sup
  if(T2 > T.sup) {
    T2 <- T.sup
    warning("Parameter 'T2' set to maximum value 'length(x) %/% 2' (",T.sup,")")
  }
  nk <- T2-T1+1
  pidx <- if(missing(nopermute)) 1L:length(x) else (1L:length(x))[-nopermute]
  #
  a <- system.time({
  	res <- .C("WRperiodogram",as.double(x),n,as.integer(T1),as.integer(T2),double(nk),
            as.integer(nperm),pidx-1L,length(pidx),integer(nk),NAOK = TRUE)[c(5L,9L)]
	})
  a[3] <- sprintf("%2f",a[3])
  if(print.time) cat("Time to compute periodogram =",a[3]," sec",'\n')
  #
  out <- matrix(NA,nk,5)
  colnames(out) <-c("Period","WR.stat","p.value","p.corrected","p.corr.progr")
  out[,1L] <- T1:T2
  out[,2L] <- sqrt(res[[1L]])
  out[,3L] <- (res[[2L]]+1)/(nperm+1)
  if(mult=="sidak") {
    out[,"p.corrected"] <- 1 - (1 - out[,"p.value"])^nk
    out[,"p.corr.progr"] <- 1 - (1 - out[,"p.value"])^(1L:nk)
  } else {
    out[,"p.corrected"] <- out[,"p.value"]*nk
    out[,"p.corr.progr"] <- out[,"p.value"]*(1L:nk)
    out[out[,"p.corrected"]>1,"p.corrected"] <- 1
    out[out[,"p.corr.progr"]>1,"p.corr.progr"] <- 1
  }
  return(structure(out,class="WRperio"))
}
#
plot.WRperio <- function(x, prog=1, alpha=0.05, line.col="red", ...)
{
  plot(x[,1], x[,2], col="red", xlab="Period", 
       ylab="Whittaker-Robinson stat.", ylim=c(0, max(x[,2])))
  lines(x[,1], x[,2], col=line.col)              # Plot W-R statistics
  points(x[,1], x[,2], pch=22, bg="white")       # Re-plot points as white squares
  if(is.na(x[1,3])) cat("p-values in the W-R periodogram were not computed\n")
  select <- which(x[,(prog+2)] <= alpha)
  points(x[select,1], x[select,2], pch=22, bg="black") # Significant points: black
}
#
