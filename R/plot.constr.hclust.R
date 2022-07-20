## **************************************************************************
##
##    (c) 2018-2022 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    **Plotting method for the constr.hclust-class**
##
##    This file is part of constr.hclust
##
##    constr.hclust is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    constr.hclust is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with constr.hclust. If not, see <https://www.gnu.org/licenses/>.
##
## /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
## |                                                            |
## |  CONSTRAINED HIERARCHICAL CLUSTERING                       |
## |  using (user-specified) criterion                          |
## |                                                            |
## |  C implementation of the Lance and Williams (1967) method  |
## |  for hierarchical clustering with or without the spatial   |
## |  contiguity constraint.                                    |
## |                                                            |
## \-----------------------------------------------------------*/
##
##    R source code file
##
## **************************************************************************
##
#' Plotting Method For Space- And Time-Constrained Clustering
#' 
#' Method \code{plot.constr.hclust} displays the results of space-constrained or
#' time-constrained agglomerative cluster analyses obtained from multivariate
#' dissimilarity matrices.
#' 
#' @usage \method{plot}{constr.hclust}(x, k, xlim, ylim, xlab, ylab, links,
#' points=TRUE, pch=21L, hybrids=c("change","single","none"), lty.hyb=1L,
#' lwd.hyb=1, col.hyb="black", plot=TRUE, col, axes, cex=1, lty, lwd, lwd.pt=1,
#' invert.axes=FALSE, ...)
#' 
#' @param x A \code{\link{constr.hclust-class}} object.
#' @param k The number of clusters to delineate.
#' @param xlim Limits, in abscissa, of the zone to be plotted.
#' @param ylim Limits, in ordinate, of the zone to be plotted.
#' @param xlab Labels for x axis annotation.
#' @param ylab Labels for y axis annotation.
#' @param links Should segments be drawn to represent the edges (links)
#' (default: FALSE).
#' @param points Should observation points be drawn (default: TRUE).
#' @param pch Point character to display observations (default: 21, a circle
#' with a background color).
#' @param hybrids How should hybrid segments be drawn (default: "change").
#' @param lty.hyb Line type to use for hybrid segments (default: lty).
#' @param lwd.hyb Width of hybrid segments with respect to lwd (default: 1).
#' @param col.hyb Colour of hybrid segments, when applicable (default: "black").
#' @param plot Should a new plotting window be opened first (default: TRUE).
#' @param col Colours to use for the \code{k} different clusters (see details). 
#' Default: for \code{k <= 10}: a set of 10 contrasted colours; otherwise, a set
#' of rainbow colours).
#' @param cex Text and symbol magnification (see \link{graphical parameters})
#' (default: 1).
#' @param axes Should the axes be displayed (default: TRUE).
#' @param lty Reference line type (see \link{graphical parameters} for details).
#' @param lwd Reference line width (see \link{graphical parameters} for
#' details).
#' @param lwd.pt Line width around points with respect to lwd (default: 1).
#' @param invert.axes Should axes be inverted on the plot (default: FALSE).
#' @param ... Other \link{graphical parameters}.
#' 
#' @details The plotting method uses the coordinates provided by the user of
#' \code{\link{constr.hclust}} to display the observations. It cuts the tree
#' (see \link{cutree}) into \code{k} clusters and uses the colours provided by
#' the user as argument \code{col} to display each cluster using the indices
#' returned by \code{\link{cutree}}. When \code{links = TRUE}, each edge is
#' displayed as a segments with colours corresponding to the clusters at its two
#' ends. A special treatment is done for hybrids edges: those whose ends lie in
#' different clusters; it is controlled by argument \code{hybrids}. When
#' argument \code{hybrids="change"} (the default), hybrid links are represented
#' as segments whose colours change halfway. When \code{hybrids="single"},
#' hybrid edges are shown as single-color lines, whose color is given as
#' argument \code{col.hyb}, whereas \code{hybrids="none"} suppresses the drawing
#' of hybrid edges. Whenever hybrid edges are displayed, their width with
#' respect to the lwd value is controlled by argument \code{lwd.hyb}.
#' 
#' When argument \code{plot=FALSE}, no \code{plot} command is issued and the
#' points (and segments when \code{links = TRUE}) are drawn over an existing
#' plotting window. This functionality is to allow one to plot the result of a
#' constrained clustering over an existing map. In that case, arguments
#' \code{xlim}, \code{ylim}, \code{axes}, and all other
#' \link{graphical parameters} to which the method \link{plot} would responds
#' are ignored.
#' 
#' The default colours are generated by function \link{rainbow}; see
#' \link{palette} for further details on using colour palettes in R. 
#' The colour palette can be changed by the user.
#' 
#' When disjoint clusters are present (i.e., when the graph provided to
#' \code{\link{constr.hclust}} is not entirely connected), the function does not
#' allow one to plot fewer clusters than the number of disjoint subsets; a
#' warning message is issued to notify the user.
#' 
#' @author Guillaume Guénard \email{guillaume.guenard@umontreal.ca}
#' and Pierre Legendre \email{pierre.legendre@@umontreal.ca}
#' 
#' @examples
#' 
#' ## Artificial map data from Legendre & Legendre (2012, Fig. 13.26)
#' ## n = 16
#' 
#' dat <- c(41,42,25,38,50,30,41,43,43,41,30,50,38,25,42,41)
#' coord.dat <- matrix(c(1,3,5,7,2,4,6,8,1,3,5,7,2,4,6,8,
#'                       4.4,4.4,4.4,4.4,3.3,3.3,3.3,3.3,
#'                       2.2,2.2,2.2,2.2,1.1,1.1,1.1,1.1),16,2)
#' 
#' ## Obtaining a list of neighbours:
#' library(spdep)
#' listW <- nb2listw(tri2nb(coord.dat), style="B")
#' links.mat.dat <- listw2mat(listW)
#' neighbors <- listw2sn(listW)[,1:2]
#' 
#' ## Calculating the (Euclidean) distance between points:
#' D.dat <- dist(dat)
#' 
#' ## Display the points:
#' plot(coord.dat, type='n',asp=1)
#' title("Delaunay triangulation")
#' text(coord.dat, labels=as.character(as.matrix(dat)), pos=3)
#' for(i in 1:nrow(neighbors))
#'     lines(rbind(coord.dat[neighbors[i,1],],
#'           coord.dat[neighbors[i,2],]))
#' 
#' ## Clustering with a contiguity constraint described by a list of
#' ## links:
#' grpWD2cst_constr_hclust <-
#'     constr.hclust(
#'         D.dat, method="ward.D2",
#'         neighbors, coord.dat)
#' 
#' ## Plot the results with k=5 clusters on a map:
#' plot(grpWD2cst_constr_hclust, k=5, links=TRUE, las=1,
#'      xlab="Eastings", ylab="Northings", cex=3, lwd=3)
#' 
#' ## Repeat the plot with other values of k (number of groups)
#' 
#' @importFrom grDevices dev.cur
#' @importFrom graphics par
#' 
#' @evalNamespace "S3method(plot,constr.hclust)" ## Waiting for a better way...
plot.constr.hclust <- function(x, k, xlim, ylim, xlab, ylab, links=FALSE,
                               points=TRUE, pch=21L,
                               hybrids=c("change","single","none"), lty.hyb=1L,
                               lwd.hyb=1, col.hyb="black", plot=TRUE,
                               col, axes=TRUE, cex=1, lty, lwd, lwd.pt=1,
                               invert.axes=FALSE, ...) {
    hybrids <- match.arg(hybrids)
    if(missing(lty)) lty <- par()$lty
    if(missing(lwd)) lwd <- par()$lwd
    
    if(!plot&&(dev.cur()==1L))
        stop("Use 'plot=FALSE' only for drawing over as existing plot!")
    if(is.null(x$coords)) {
        class(x) <- "hclust"
        plot(x, cex=cex, ...)
        return(invisible(NULL))
    }
    if(missing(col))
      if(k > 10) {
        col <- rainbow(1.2*k)[1L:k]
      } else
        col <- c("blue", "gold", "grey70", "cadetblue2", "red", "orange3",
                 "coral2", "green", "blueviolet", "grey30")[1L:k]
    if(any(nna <- is.na(x$height)))
        if(k < sum(nna) + 1L) {
            warning("Impossible to plot the cluster for k = ", k, " because ",
                    "the graph provided involves ", sum(nna) + 1L, " non",
                    "-connected clusters that cannot be linked at any ",
                    "dissimilarity level. This cluster is plotted for k = ",
                    sum(nna) + 1L)
            k <- sum(nna) + 1L
        }
    cl <- cutree(x,k)
    coords <- x$coords
    if(invert.axes) coords <- coords[,2L:1L]
    if(plot) {
        if(missing(xlim)) {
            xlim <- range(coords[,1L])
            if((xlim[2L]-xlim[1L])==0)
                xlim <- xlim + c(-0.5,0.5)
        }
        if(missing(ylim)) {
            ylim <- range(coords[,2L])
            if((ylim[2L]-ylim[1L])==0)
                ylim <- ylim + c(-0.5,0.5)
        }
        if(missing(xlab))
            xlab <- if(diff(range(coords[,1L]))==0) "" else "x"
        if(missing(ylab))
            ylab <- if(diff(range(coords[,2L]))==0) "" else "y"
        plot(NA, asp=1, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, axes=axes,
             type="n", cex=cex, ...)
    }
    if(links) {
        if(is.null(x$links))
            x$links <- cbind(1L:(nrow(coords)-1L),2L:nrow(coords))
        dev.hold()
        for(i in 1:nrow(x$links)) {
            ## i=1L
            clij <- cl[x$links[i,1L:2L]]
            if(clij[1L]==clij[2L]) {
                segments(coords[x$links[i,1L],1L],
                         coords[x$links[i,1L],2L],
                         coords[x$links[i,2L],1L],
                         coords[x$links[i,2L],2L],
                         col=col[clij[1L]], lty=lty, lwd=lwd, ...)
            } else {
                ## The link is an hybrid
                if(hybrids=="change") {
                    mid <- c(mean(coords[x$links[i,],1L]),
                             mean(coords[x$links[i,],2L]))
                    segments(coords[x$links[i,1L],1L],
                             coords[x$links[i,1L],2L],
                             mid[1L], mid[2L],
                             col=col[clij[1L]],
                             lty=lty.hyb,
                             lwd=lwd*lwd.hyb, ...)
                    segments(mid[1L], mid[2L],
                             coords[x$links[i,2L],1L],
                             coords[x$links[i,2L],2L],
                             col=col[clij[2L]],
                             lty=lty.hyb,
                             lwd=lwd*lwd.hyb, ...)
                } else if(hybrids=="single")
                    segments(coords[x$links[i,1L],1L],
                             coords[x$links[i,1L],2L],
                             coords[x$links[i,2L],1L],
                             coords[x$links[i,2L],2L],
                             col=col.hyb,
                             lty=lty.hyb,
                             lwd=lwd*lwd.hyb, ...)
            }
        }
        dev.flush()
    }
    if(points)
        points(x=coords[,1L], y=coords[,2L], bg=col[cl], lty=lty,
               lwd=lwd.pt*lwd, cex=cex, pch=pch, ...)
    return(invisible(NULL))
}
## 
