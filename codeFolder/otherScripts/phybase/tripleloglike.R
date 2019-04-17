tripleloglike<-function(sptree,spname,dna)
{
	ntaxa<-dim(dna)[1]
	ntriple<-ntaxa*(ntaxa-1)*(ntaxa-2)/6
	triple<-matrix(0,nrow=ntriple,ncol=5)
	par<-matrix(0,nrow=ntriple,ncol=4)
	nodematrix<-read.tree.nodes(sptree,name=spname)$node

	triple<-triplenumber(dna)

	for(i in (ntaxa+1):(2*ntaxa-2))
	{
		son1<-nodematrix[i,2]
		sonnode1<-offspring.species(son1, nodematrix, ntaxa)
		son2<-nodematrix[i,3]
		sonnode2<-offspring.species(son2, nodematrix, ntaxa)
		father <- i

		while(father != (2*ntaxa-1))
		{
			son <- father
			father <- nodematrix[father,1]
			{if(nodematrix[father,2] == son)
				sonnode3 <- offspring.species(nodematrix[father,3],nodematrix,ntaxa)
			else
				sonnode3 <- offspring.species(nodematrix[father,2],nodematrix,ntaxa)}

			

			for(j in 1:length(sonnode1))
				for(k in 1:length(sonnode2))
					for( w in 1:length(sonnode3))
					{	
						par[sonnode1[j]+sonnode2[k]+sonnode3[w]-5,] <- triplepara(i,father,nodematrix,ntaxa)	
					}
		}
	}

	triplep<-apply(par,1,tripleProb) 
	loglike <- sum(triple*log(t(triplep)))

	loglike

}
