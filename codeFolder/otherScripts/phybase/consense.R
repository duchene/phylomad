`consense` <-
function(treestr,name,type="freq")
{

    if(length(grep("freq",type,ignore.case=TRUE))==0 & length(grep("prop",type,ignore.case=TRUE))==0)
	stop("type is either freq or prop!")
    if(length(grep("freq",type,ignore.case=TRUE))==1)
	freq<-1  else  freq<-0
      

   

    ntree<-length(treestr)
    tree<-read.tree.nodes(treestr[1],name)$nodes
    nspecies<-length(name)  ##rooted tree
    con<-partition.tree(tree,nspecies)
    

    ##partition trees
    for(i in 2:ntree){
	tree<-read.tree.nodes(treestr[i],name)$nodes
	partition<-partition.tree(tree,nspecies)

 	##frequency of the partition
	for(j in 1:dim(partition)[1]){
		found<-FALSE
		for(k in 1:dim(con)[1]){
		 	if(sum(abs(partition[j,1:nspecies]-con[k,1:nspecies]))==0){
				con[k,nspecies+1]<-con[k,nspecies+1]+1
				found<-TRUE
				break
			} 
		}
		if(!found)
			con<-rbind(con,partition[j,])
	}
		
   	
    }
    con<-con[order(con[,(dim(con)[2])],decreasing=TRUE),]
    majcon<-con[con[,nspecies+1]>ntree/2,]
    majcon<-majcon[order(apply(majcon[,1:nspecies],1,sum)),]
    string<-paste(1:nspecies)
    round<-dim(majcon)[1]

    if(sum(majcon[round,1:nspecies]) != nspecies){
	newrow<-rep(1,nspecies+1)
	newrow[nspecies+1]<-ntree
	majcon<-rbind(majcon,newrow)
        round<-round+1
    }

    
    for(i in 1:round){

        ncol<-dim(majcon)[2]-1
  	b<-majcon[i,1:ncol]*(1:ncol)
	b<-b[b>0]
   	
	##update groups
	newname<-paste(string[b],sep="",collapse=",")
    	newname<-paste("(",newname,sep="")
	newname<-paste(newname,")",sep="")
	if(i != round){	
    		newname<-paste(newname,":",sep="")
		if(freq)
			newname<-paste(newname,majcon[i,ncol+1],sep="")
		else
			newname<-paste(newname,round(majcon[i,ncol+1]/ntree,2),sep="")
	}


    	string[b[1]]<-newname
	
	##update dist,string, and branch length
	index<-1:(ncol+1)
	index[b]<-0
	index[b[1]]<-b[1]
	index<-index[index>0]
        
    	majcon<-majcon[1:dim(majcon)[1],index]
        string<-string[index[1:(length(index)-1)]]
    }


    z<-list(treestr=as.character,name=as.character)
    z$treestr<-string
    z$name<-name 
    z  
}