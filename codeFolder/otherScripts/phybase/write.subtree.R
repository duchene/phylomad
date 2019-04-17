`write.subtree` <-
function (inode, nodematrix, taxaname, root)
{
	if(nodematrix[inode,2] > 0)
	{
		son1 <- nodematrix[inode,2]
		son2 <- nodematrix[inode,3]
		x <- paste("(",write.subtree(son1,nodematrix,taxaname, root),":",nodematrix[son1,4],",",write.subtree(son2,nodematrix,taxaname,root),":",nodematrix[son2,4],")",sep="")		
	}
	else x <- taxaname[inode]
	if(inode == root) x <- paste(x,";",sep="")
	return(x)
}

