'getcoaltime'<-
function(genetree,sptree,ntax,nspecies,species.structure)
{
	geneanc<-rep(0,dim(genetree)[1])
	index<-1
	coal<-matrix(-1,nrow=(dim(genetree)[1]-ntax),ncol=2)
	for(i in 1:nspecies)
	{
		genenode<-which(species.structure[i,]==1)
		spancs<-ancandtime(i,sptree,nspecies)
		
		for(j in 1:length(genenode))
		{
			inodegene<-genenode[j]
			geneanc1<-ancandtime(inodegene,genetree,ntax)
	
			geneanc2<-rep(0,dim(genetree)[1])
			geneanc2[geneanc1[,1]]<-1
			newanc<-1:sum((geneanc2-geneanc)>0)
		
			coaltime<-geneanc1[newanc,2]
			for(k in 1:length(coaltime))
			{
			if(k != 1){
			a<-which(spancs[,2]<=coaltime[k])
			coal[index,1]<-spancs[length(a),1]
			coal[index,2]<-coaltime[k]-spancs[a[length(a)],2]
			index<-index+1}
			}
			geneanc[geneanc1[newanc,1]]<-1
		}

	}
	return(coal)
}