mutation_exp<-function(sptree,root,inode,nspecies,alpha)
{
	if(inode == root)
	{
		sptree[root,6]<-1.0
		sptree<-mutation_exp(sptree,root,sptree[root,2],nspecies,alpha)
		sptree<-mutation_exp(sptree,root,sptree[root,3],nspecies,alpha)
	}
	else
	{
		sptree[inode,6]<-rgamma(1,alpha,(alpha/sptree[sptree[inode,1],6]))
		if(inode > nspecies)
		{
			sptree<-mutation_exp(sptree,root,sptree[inode,2],nspecies,alpha)
			sptree<-mutation_exp(sptree,root,sptree[inode,3],nspecies,alpha)
		}
	}
	return(sptree)
}
