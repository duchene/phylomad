`maxtree` <-
function(genetreevector,spname,taxaname,species.structure)
{
	ntree<-length(genetreevector)
    	nodes1<-read.tree.nodes(genetreevector[1],taxaname)$nodes
    	nspecies<-length(spname)
	ntaxa<-length(taxaname)

    	dist<-pair.dist(nodes1,ntaxa)
    	#treestr<-paste(1:nspecies)
    	#brlens<-rep(0,nspecies)
    
   	#calculate pairwise distance
    	for(i in 2:ntree)
	{
		nodes1<-read.tree.nodes(genetreevector[i],taxaname)$nodes
		dist1<-pair.dist(nodes1,ntaxa)
		dist<-pmin(dist,dist1)
    	}
	
	spdist1<-matrix(0,nrow=nspecies,ncol=ntaxa)
	for(i in 1:nspecies)
		{
			taxa<-species.structure[i,]*(1:ntaxa)
			taxa<-taxa[taxa>0]
			if(length(taxa)>1)  spdist1[i,]<-apply(dist[taxa,],2,min)
			else spdist1[i,]<-dist[taxa,]
		}

	spdist<-matrix(0,nrow=nspecies,ncol=nspecies)
	for(i in 1:nspecies)
		{
			taxa<-species.structure[i,]*(1:ntaxa)
			taxa<-taxa[taxa>0]
			if(length(taxa)>1)  spdist[,i]<-apply(spdist1[,taxa],1,min)
			else spdist[,i]<-spdist1[,taxa]
		}

	dist<-spdist/2

	treestr<-upgma(dist,spname,method="min")$treestr

	
	return(treestr)
	
}