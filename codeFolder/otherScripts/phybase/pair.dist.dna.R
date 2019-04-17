`pair.dist.dna` <-
function(sequences,nst=0)
{
    
    nsequences<-dim(sequences)[1]
    seqlength<-dim(sequences)[2]

    dist<-matrix(0,nsequences,nsequences)
    for(i in 1:(nsequences-1)){
	for(j in (i+1):nsequences){
		dist[i,j]<-(seqlength-sum(sequences[i,]==sequences[j,]))/seqlength
	}
    }
    if(nst==1)
	dist<--0.75*log(1-4*dist/3)
    
    return(t(dist)+dist)
}