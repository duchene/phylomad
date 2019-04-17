`swap.nodes` <-
function(inode, jnode, name, nodematrix)
{
    	nodes<-nodematrix
    	##swap father
    	father1<-nodes[inode,1]
    	father2<-nodes[jnode,1]
    	nodes[inode,1]<-father2
    	nodes[jnode,1]<-father1

    	##swap son
    	if(nodes[father1,2] == inode)
		nodes[father1,2]<-jnode 
	else 
		nodes[father1,3]<-jnode
    	if(nodes[father2,2] == jnode)
  		nodes[father2,2]<-inode 
	else 
		nodes[father2,3]<-inode 
    	z <- list(nodes = matrix(0, dim(nodes)[1], dim(nodes)[2]), treestr="")
    	z$nodes<-nodes
    	z$treestr<-subtree(dim(nodes)[1],name,nodes)
    	z 
}

