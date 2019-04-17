`mrca.2nodes` <-
function(inode,jnode,nodematrix)
{
    ancestor1<-ancestor(inode,nodematrix)
    ancestor2<-ancestor(jnode,nodematrix)

    for(i in 1:length(ancestor1)){
	s<-FALSE
	for(j in 1:length(ancestor2)){
		if((ancestor2[j]-ancestor1[i])==0){
			s<-TRUE
  			break
		}
    	}
 	if(s)
		break
    }
    mrca<-ancestor1[i]

    totallength<-0
    if(i > 1){	
    leftlength<-nodematrix[ancestor1[1:(i-1)],4]
    leftlength<-leftlength[leftlength>0]
    totallength<-sum(leftlength)
    } 
    if(j>1)
    {
    rightlength<-nodematrix[ancestor2[1:(j-1)],4]
    rightlength<-rightlength[rightlength>0]
    totallength<-totallength+sum(rightlength)
    }
    z<-list(anc=as.integer,dist=as.double)
    z$anc<-mrca
    z$dist<-totallength
     
    z
  
}
