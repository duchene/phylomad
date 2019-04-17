`subtree` <-
function(inode, name, nodematrix)
{
    	if(inode<1 | inode>dim(nodematrix)[1])
	{
		a<-paste("The node number should be between 1 and",dim(nodematrix)[1])
		stop(a)
    	}
    	nspecies<-length(name) 
  
    	if(inode<=nspecies)
	{
       	treestr<-paste("The tree has only one taxon: ",inode)
		return(treestr)
    	}
    	treestr<-write.subtree(inode, nodematrix, nspecies,inode)
    	z<-node2name(treestr,name)
    	z  
}

