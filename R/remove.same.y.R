`remove.same.y` <-
function(coords,link){
	
	#### Remove links with the same y coordinates (after rotation)
	#### Find all the points (sites) with the same y coordinate
	#CC# Make sure the names match
	colnames(link)<-c("from","to")
	
	#CC# Creat object link.tmp and coords.tmp
	link.tmp<-link
	coords.tmp<-coords
	
	#CC# Loop checking for sites with the same y coordinates
	while(nrow(coords.tmp)!=0){
		#CC# Which site has the same coordinates
		same.coords<-which(coords.tmp[,3]==coords.tmp[1,3])
		#CC# Find the sites numbers which have the same y coordinates
		same.y<-coords.tmp[same.coords,1]
		#CC# What to do if there are more than one site with the same y coordinates
		if(length(same.y) > 1){
			#CC# Construct a matrix in which will be entered links with the same
			#CC# y coordinates as i in the "from" section.
			same.from<-matrix(c(0,0),nrow=1,ncol=2)
			#CC# Name each coloumn of the site for smoother running of the function
			colnames(same.from)<-c("from","to")
			#CC# Loop adding a row for each link with the same y coordinates in the 
			#CC# "from" section
			for(j in same.y){
				#CC# Which sites have the same "from"
				same.from<-rbind(same.from,as.data.frame(link[which(link[,1]==j),]))
			}
			#CC# Remove the 0 0 which was there only to start the matrix
			same.from<-same.from[-1,]
			#CC# Construct a matrix in which will be entered links with the same
			#CC# y coordinates as i in the "to" section
			same<-matrix(c(0,0),nrow=1,ncol=2)
			#CC# Name each coloumn of the site for smoother running of the function
			colnames(same)<-c("from","to")
			#CC# Loop adding a row for each link with the same y coordinates in the 
			#CC# "to" section
			for(j in same.y){
				#CC# Which sites have the same "to"
				same<-rbind(same,as.data.frame(same.from[which(same.from[,2]==j),]))
			}
			#CC# Remove the 0 0 which was there only to start the matrix
			same<-same[-1,]
			#CC# If the object "same" has a positive number of rows
			if(nrow(same)!=0){
				#### Put all the site in a single line (this will limit the repetition)
				same.no.rep<-c(same[,1],same[1:nrow(same),2])
				
				#### To remove the horizontal links
				#CC# Match the links in "same" with the links in "link" in the first
				#CC# column (from)
				link.from<-which(link.tmp[,1] %in% same.no.rep)
				#CC# Match the links in "same" with the links in "link" in the second
				#CC# column (to)
				link.to<-which(link.tmp[,2] %in% same.no.rep)
				#CC# See which from link.from and link.to have the same position
				link.from.to<-link.from[which(link.from%in%link.to)]
				#CC# Remove the link that have the same position in "link"
				link.tmp<-link.tmp[-link.from.to,]
			
				#CC# Remove the coordinates of the site considered
				coords.tmp<-as.data.frame(coords.tmp[-same.coords,])
			}else{
				#CC# Remove the coordinates of the site considered
				coords.tmp<-as.data.frame(coords.tmp[-same.coords,])
			}
		}else{
			#CC# Remove the coordinates of the site considered
			coords.tmp<-as.data.frame(coords.tmp[-same.coords,])
		}
	}
	#CC# Rename the object "link.tmp" by "link"
	return(link.tmp)
	
}

