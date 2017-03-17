#' Rotate a set of point by a certain angle
#' 
#' Rotate a set of XY coordinates by an angle (in radians)
#' 
#' @param xy A 2-columns matrix or data frame containing a set of X and Y
#'   coordinates.
#' @param angle Numeric. A scalar giving the angle at which the points should be
#'   rotated. The angle is in radians.
#'   
#' @return
#' 
#' A 2-columns matrix of the same size as \code{xy} giving the rotated
#' coordinates.
#' 
#' @author F. Guillaume Blanchet
#'   
#' @examples
#' 
#' ### Create a set of coordinates
#' coords<-cbind(runif(20),runif(20))
#' 
#' ### Create a series of angles
#' rad<-seq(0,pi,l=20)
#' 
#' for(i in rad){
#' 	coords.rot<-rotation(coords,i)
#' 	plot(coords.rot)
#' }
#' 
#' ### Rotate the coordinates by an angle of 90 degrees
#' coords.90<-rotation(coords,90*pi/180)
#' coords.90
#' 
#' plot(coords,xlim=range(rbind(coords.90,coords)[,1]),ylim=range(rbind(coords.90,coords)[,2]),asp=1)
#' points(coords.90,pch=19)
#' 
#' @keywords manip
#' @export

`rotation` <-
    function(xy, angle)
    {
        xy <- as.matrix(xy)
        ### Find cos and sin of the angle
        cos.angle <- cos(angle)
        sin.angle <- sin(angle)
        
        ### Rotate the set of coordinates
        xy.rot <- xy %*% t( matrix(c(cos.angle,sin.angle, -sin.angle, cos.angle), 2,2) )
        
        return(xy.rot)
    }

