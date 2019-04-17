"ChangeBrlen"<-
function(sptree, spnodedepth, spnode, genetree, genenodedepth, genenode)
{
	inode <- FindSpnodeDownGenenode(sptree,spnodedepth,spnode,genenodedepth,genenode)
	jnode <- FindSpnodeDownGenenode(sptree,spnodedepth,spnode,genenodedepth,genetree[genenode,1])

	if(inode == jnode)
	{
		length <- (genetree[genenode,4]) * sptree[inode,6]
	}
	else
	{
		father <- sptree[inode,1]
		length <- (spnodedepth[father] - genenodedepth[genenode])*sptree[inode,6]
		while(father != jnode)
		{
			inode <- father;
			father <- sptree[father,1] 
			length <- length + (spnodedepth[father] - spnodedepth[inode])*(sptree[inode,6])
		}
		length <- length + (genenodedepth[genetree[genenode,1]] - spnodedepth[father])*sptree[father,6]
  	}
	return(length)		
}