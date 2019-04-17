`ancestor` <-
function(inode,nodematrix)
{
	if(!is.rootedtree(nodematrix))
	{
		warnings("The tree is not rooted!")
	}
    	rootnode<-rootoftree(nodematrix)
	nnodes<-dim(nodematrix)[1]

    	if(inode == rootnode)
		ancestor<-rootnode 
	else
	{
    		ancestor<-rep(0,nnodes) 
    		ancestor[1]<-inode    
    		i<-1
    		while(ancestor[i] != rootnode)
		{
       	ancestor[i+1]<-nodematrix[ancestor[i],1]
		i<-i+1
    		}
    }
    return(ancestor[ancestor>0])
}

