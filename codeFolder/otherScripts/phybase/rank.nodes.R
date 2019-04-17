`rank.nodes` <-
function(treenode, inode, ntaxa, start, rank)
{   
    	if(inode > ntaxa)
    	{
		left<-treenode[inode,2]
		right<-treenode[inode,3]
		rank[inode]<-start
		rank[left]<-start-1
		rank[right]<-start-1
		rank<-rank.nodes(treenode,left,ntaxa,rank[left],rank)
		rank<-rank.nodes(treenode,right,ntaxa,rank[right],rank)
		return(rank)
    	}
    	else
    	{
		return(rank)
    	} 
}
