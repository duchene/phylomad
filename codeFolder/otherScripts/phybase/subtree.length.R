`subtree.length` <-
function(inode,nodes,nspecies)
{
    	if(inode<=nspecies)
		length<-0
    	if(inode>nspecies)
	{
		son1<-nodes[inode,2]
        	son2<-nodes[inode,3]
		length<-subtree.length(son1,nodes,nspecies)+subtree.length(son2,nodes,nspecies)+nodes[son1,4]+nodes[son2,4]
    	}
    	return(length)
}

