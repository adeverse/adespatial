`rm.double.link` <-
function (link) 
{
    for (i in 1:nrow(link)) {
        tmp <- which(link[i, 1] == link[, 2] & link[i, 2] == 
            link[, 1])
        if (length(tmp) != 0) {
            link[tmp, ] <- c(0, 0)
        }
    }
    tmp <- which(link[, 1] != 0 & link[, 2] != 0)
    if (length(tmp) != 0) 
        link <- link[tmp, ]
    return(link)
}

