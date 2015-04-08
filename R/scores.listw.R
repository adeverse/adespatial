scores.listw <- function (listw,  wt = rep(1, length(listw$neighbours)), MEM.autocor = c("non-null", "all", "positive", 
    "negative")) 
{
    if (!inherits(listw, "listw")) 
        stop("not a listw object")
    MEM.autocor <- match.arg(MEM.autocor)

    wt <- wt / sum(wt)
    w <- listw2mat(listw)
    sumW <- sum(w)
    n <- nrow(w)
    symmetric <- isSymmetric.matrix(w, check.attributes = FALSE)
    if (symmetric == FALSE) {
        w <- (w + t(w))/2
    }
    w <- bicenter.wt(w, row.wt = wt, col.wt = wt)
    res <- eigen(w, symmetric = TRUE)
    eq0 <- which(apply(as.matrix(res$values/max(abs(res$values))), 
        1, function(x) identical(all.equal(x, 0), TRUE)))
    if (length(eq0) == 0) {
        stop(" Illegal matrix: no null eigenvalue")
    }
   
    if(MEM.autocor == "all"){
        if (length(eq0) == 1) {
            res$values <- res$values[-eq0]
            res$vectors <- res$vectors[, -eq0]
            
        } else if (length(eq0) > 1) {
            w <- cbind(rep(1, n), res$vectors[, eq0])
            w <- qr.Q(qr(w))
            res$values[eq0] <- 0
            res$vectors[, eq0] <- w[, -ncol(w)]
            res$values <- res$values[-eq0[1]]
            res$vectors <- res$vectors[, -eq0[1]]
        }
    } else if(MEM.autocor == "non-null"){        
        res$values <- res$values[-eq0]
        res$vectors <- res$vectors[, -eq0]
    } else if (MEM.autocor == "positive") {
        posi <- which(res$values > -sumW/((n - 1) * n))
        res$values <- res$values[posi]
        res$vectors <- res$vectors[, posi]
    } else if (MEM.autocor == "negative") {
        neg <- sort(which(res$values < -sumW/((n - 1) * n)), 
                    decreasing = TRUE)
        res$values <- res$values[neg]
        res$vectors <- res$vectors[, neg]
    }

    a0 <- as.data.frame(res$vectors) / sqrt(wt)
    z <- res$values
    row.names(a0) <- attr(listw,"region.id")
    names(a0) <- paste("MEM", 1:ncol(a0), sep = "")
    attr(a0,"values") <- z
    attr(a0,"weights") <- wt
    attr(a0,"call") <- match.call()
    attr(a0,"class") <- c("orthobasisSp", "orthobasis","data.frame")
 
    return(a0)
}

mem <- orthobasis.listw <- scores.listw
