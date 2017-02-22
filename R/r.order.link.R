`r.order.link` <-
function (nrow.link, link, coords) 
{
    for (i in 1:nrow.link) {
        if (link[i, 2] != 0) {
            if (coords[link[i, 3], 3] < coords[link[i, 2], 3]) {
                link[i, 2:3] <- c(link[i, 3], link[i, 2])
            }
        }
    }
    return(link)
}

