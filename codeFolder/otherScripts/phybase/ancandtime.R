'ancandtime'<-
function (inode, nodematrix,nspecies) 
{
    if (!is.rootedtree(nodematrix)) {
        warnings("The tree is not rooted!")
    }
    rootnode <- rootoftree(nodematrix)
    
    nnodes <- dim(nodematrix)[1]
    if (inode == rootnode) 
        ancestor <- rootnode
    else {
        ancestor <- matrix(0, nrow=nnodes,ncol=2)
        ancestor[1,1] <- inode
	ancestor[1,2] <- node.height(inode,nodematrix,nspecies)
        i <- 1
        while (ancestor[i] != rootnode) {
            ancestor[i+1,1] <- nodematrix[ancestor[i,1], 1]
	    ancestor[i+1,2] <- ancestor[i,2]+nodematrix[ancestor[i,1],4]
            i <- i + 1
        }
    }
    anc<-ancestor[which(ancestor[,1]>0),]
    colnames(anc)<-c("anc","time")
    anc
}