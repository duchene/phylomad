"sim.coaltree.sp.mu"<-
function(sptree, spname,seq,numgenetree,method="dirichlet",alpha=5.0)
{
	nodematrix<-read.tree.nodes(sptree,spname)$nodes
	rootnode<-dim(nodematrix)[1]
	nspecies<-(rootnode+1)/2
	ntaxa<-sum(seq)
	
	#generate mutation rates
	#nodematrix<-cbind(nodematrix,rep(-100,rootnode))
	if(tolower(method) == "gamma")
		nodematrix<-mutation_exp(nodematrix,rootnode,rootnode,nspecies,alpha)
	if(tolower(method) == "dirichlet")
		nodematrix[,6] <- rdirichlet(1,rep(alpha,dim(nodematrix)[1]))*dim(nodematrix)[1]
	if(tolower(method) == "user")
		nodematrix[,6] <- alpha

	index<-1
	seqname<-rep("",ntaxa)
	for(i in 1:nspecies)
		for(j in 1:seq[i])
		{
			if(seq[i] > 1)
				seqname[index]<-paste(spname[i],"s",j,sep="")
			else
				seqname[index]<-spname[i]
			index<-index+1
		}

	speciesmatrix<-matrix(0,nrow=nspecies,ncol=ntaxa)

	index<-1	
	for(i in 1:length(seq))
	{
		for(j in 1:seq[i])
		{
			speciesmatrix[i,index]<-1
			index<-index+1
		}
	}
	
	spnodedepth<-rep(0,2*nspecies-1)
	for(i in 1:(2*nspecies-1))
	{
		spnodedepth[i]<-node.height(i,nodematrix,nspecies)
	}

	treestr<-rep("",numgenetree)
	for(j in 1:numgenetree)
	{
		str<-sim.coaltree.sp(rootnode,nodematrix,nspecies,seq,name=spname)$gt
		genetree<-read.tree.nodes(str,name=seqname)$nodes
		genenodedepth<-rep(0,2*ntaxa-1)
	
		for(i in 1:(2*ntaxa-1))
		{
			genenodedepth[i]<-node.height(i,genetree,ntaxa)
		}	
		coaltree<-populationMutation(nodematrix,spnodedepth,genetree,genenodedepth,speciesmatrix)
		treestr[j]<-write.subtree(dim(coaltree)[1],coaltree,seqname,dim(coaltree)[1])
	}
	z <- list(gt=as.character, st=as.matrix,seqname=as.character)
    	z$gt <- treestr
    	z$st <- nodematrix
	z$seqname<-seqname
	return(z)
}
