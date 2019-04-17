triplepara<-function(inode,jnode,nodematrix,nspecies)
{
	par<-rep(0,4)
	height1<-node.height(inode, nodematrix, nspecies)
	height2<-node.height(jnode, nodematrix, nspecies)

	if(height1 < height2)
	{
		par[1] <- height2-height1
		par[2] <- height1
		par[3] <- nodematrix[jnode,5]
		par[4] <- nodematrix[inode,5]
	}
	else if(height1 > height2)
	{
		par[1] <- height1-height2
		par[2] <- height2
		par[3] <- nodematrix[inode,5]
		par[4] <- nodematrix[jnode,5]
	}
	else
	{
		warnings("something is wrong in triplepara")
	}
	par
}