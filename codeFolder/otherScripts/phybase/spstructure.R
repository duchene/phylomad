'spstructure'<-
function(numsgenenodes)
{

	nspecies<-length(numsgenenodes)
	ntaxa<-sum(numsgenenodes)
	sp<-matrix(0,nrow=nspecies,ncol=ntaxa)
	index<-matrix(0,nrow=ntaxa,ncol=2)
	index[,2]<-1:ntaxa

	b<-1
	for(i in 1:nspecies)
	{
		for(j in 1:numsgenenodes[i])
		{
			index[b,1]<-i
			b<-b+1
		}
	}
	sp[index]<-1
	return(sp)
}