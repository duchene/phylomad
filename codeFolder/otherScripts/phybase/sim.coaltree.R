`sim.coaltree` <-
function(nspecies, theta)
{
  
    treestr<-paste(1:nspecies)
    height<-0
    brlens<-rep(0,nspecies)
    for(i in nspecies:2){

	##randomly choose two nodes
 	b<-sample(1:i,2)
	height<-height+rexp(1,rate=i*(i-1)/theta)

	##update groups
    	newname<-paste("(",treestr[b[1]],sep="")
    	newname<-paste(newname,":",sep="")
    	newname<-paste(newname,round(height-brlens[b[1]],5),sep="")
    	newname<-paste(newname,",",sep="")   
    	newname<-paste(newname,treestr[b[2]],sep="")
    	newname<-paste(newname,":",sep="")
    	newname<-paste(newname,round(height-brlens[b[2]],5),sep="")
    	newname<-paste(newname,")",sep="")
    	treestr[b[1]]<-newname
	brlens[b[1]]<-height

	##update dist,treestr, and branch length
	index<-1:i
    	index[b[2]]<-0
    	index<-index[index>0]
        treestr<-treestr[index]
	brlens<-brlens[index]
    	
    }
    treestr

}



