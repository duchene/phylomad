
'simSeqfromSp' <- 
function (sptree, spname, ntaxasp, ngene, theta=0, noclock=0, simsequence=1, murate="Dirichlet",alpha=5, seqlength=100, model=1, kappa=2, rate=c(1,1,1,1,1,1), frequency=c(1/4,1/4,1/4,1/4), outfile, format="phylip")
{
	if(!is.rootedtree(sptree))
		stop("the species tree must be rooted!")

	nspecies<-length(spname)
	ntaxa <- sum(ntaxasp)

	if(length(ntaxasp) != nspecies)
	{
		stop("The number of species in ntaxasp does not match the number of species in sptree!")
	}

	rootnode<-nspecies*2-1
	taxaname<-rep("", ntaxa)

	index <- 1
	
	for(i in 1:nspecies)
	{
		for(j in 1:ntaxasp[i])
		{
			if(ntaxasp[i] > 1)
				taxaname[index]<-paste(spname[i],"s",j,sep="")
			else
				taxaname[index]<-spname[i]
			index <- index + 1
		}
	}

	species.structure<-matrix(0,nspecies, ntaxa)
	index <- 1
	for(i in 1:nspecies)
	{
		for(j in 1:ntaxasp[i])
		{
			species.structure[i,index]<-1
			index <- index + 1
		}
	}
	

	#simulate sequences
	nodematrix<-read.tree.nodes(sptree,spname)$nodes
	if(theta > 0)
		nodematrix[,5]<-theta
	tree <- rep("", ngene)

	dnaseq <- matrix("", nrow = ntaxa, ncol = seqlength * ngene)
	partition <- matrix(0, nrow = ngene, ncol = 2)

	for(i in 1:ngene)
	{
		if(noclock == 0)
			tree[i]<-sim.coaltree.sp(rootoftree(nodematrix),nodematrix, nspecies,ntaxasp, name=spname)$gt
		else if(noclock == 1)
		{			
			tree[i]<-sim.coaltree.sp.mu(sptree, spname, ntaxasp, 1, method=murate,alpha=alpha)$gt
		}
		else
			stop ("noclock could be 0 or 1!")

		if(simsequence)
		{
			treenode<-read.tree.nodes(tree[i],taxaname)$nodes	
			treenode[2*ntaxa-1,4]<--9
			dna<-sim.dna(treenode,seqlength, model=model, kappa=kappa, rate=rate, frequency=frequency)
			dna[dna==1]<-"A"
			dna[dna==2]<-"C"
			dna[dna==3]<-"G"
			dna[dna==4]<-"T"
			dnaseq[, (1+(i-1)*seqlength):(seqlength*i)] <- dna
			if(tolower(format) == "phylip")
			{
				if(i==1) write.dna(dna,taxaname, file=outfile,format=format,append=FALSE)
   	 			else write.dna(dna,taxaname,file=outfile,format=format,append=TRUE)
			}
			partition [i, ] <- c(1+(i-1)*seqlength, seqlength*i)
		}
	}

	if(tolower(format) == "nexus")	
		write.dna(dnaseq, taxaname, file=outfile,partition = partition, format=format, taxa = ntaxasp, program="best", append=FALSE)

	z<-list(gt="")
	z$gt <- tree
	z
}