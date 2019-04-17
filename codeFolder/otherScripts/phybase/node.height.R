`node.height` <-
function(inode,nodematrix,nspecies)
{
    	if(inode<=nspecies)
		height<-0
    	if(inode>nspecies)
	{
		son1<-nodematrix[inode,2]
		height<-node.height(son1,nodematrix,nspecies)+nodematrix[son1,4]
    	}
    	return(height)
}

