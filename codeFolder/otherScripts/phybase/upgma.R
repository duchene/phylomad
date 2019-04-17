`upgma` <-
function(dist,name,method="average")
{
    nspecies<-length(name)
    treestr<-paste(1:nspecies)
    brlens<-rep(0,nspecies)
    b<-rep(0,2)

    for(i in 1:(nspecies-1)){
	##find the minimum distance
    	diag(dist)<-NA
    	position <- which(dist == min(dist, na.rm = TRUE))[1]
      b[1] <- floor(position/length(dist[1, ]))
      b[2] <- position - b[1] * length(dist[1, ])
      if (b[2] == 0) 
          	b[2] <- length(dist[1, ]) else b[1] <- b[1] + 1
	diag(dist)<-0

	##recalculate distance
	if(method == "average")
	{
    		dist[,b[1]]<-(dist[,b[1]]+dist[,b[2]])/2
    		dist[b[1],]<-(dist[b[1],]+dist[b[2],])/2
	}
	if(method == "min")
	{
		bdist<-dist[b[2],b[1]]
		dist[,b[1]]<-pmin(dist[,b[1]],dist[,b[2]])
    		dist[b[1],]<-pmin(dist[b[1],],dist[b[2],])
		dist[b[2],b[1]] <- bdist
		dist[b[1],b[2]] <- bdist
	}
 	
	##update groups
    	newname<-paste("(",treestr[b[1]],sep="")
    	newname<-paste(newname,":",sep="")
    	newname<-paste(newname,round((dist[b[2],b[1]]-brlens[b[1]]),5),sep="")
    	newname<-paste(newname,",",sep="")   
    	newname<-paste(newname,treestr[b[2]],sep="")
    	newname<-paste(newname,":",sep="")
    	newname<-paste(newname,round((dist[b[2],b[1]]-brlens[b[2]]),5),sep="")
    	newname<-paste(newname,")",sep="")
    	treestr[b[1]]<-newname
	brlens[b]<-dist[b[2],b[1]]

	##update dist,treestr, and branch length
	index <- 1:(nspecies + 1 - i)
    	index[b[2]]<-0
    	index<-index[index>0]
    	dist<-dist[index,index]
      treestr<-treestr[index]
    	brlens<-brlens[index]
    }
    treestr<-paste(treestr,";",sep="")
    z<-list(nodes=matrix,treestr="",name="")
    node<-read.tree.nodes(treestr,name)
    z$nodes<-node$nodes
    z$name<-node$name
    z$treestr<-treestr
    z
}
