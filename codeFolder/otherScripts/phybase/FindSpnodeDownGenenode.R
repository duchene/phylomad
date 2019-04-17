"FindSpnodeDownGenenode"<-
function(sptree, spnodedepth, spnode, genenodedepth, genenode)
{
	findnode<-spnode;
	root<-dim(sptree)[1]

	depth <- genenodedepth[genenode]
	father <- sptree[spnode,1]

	if(genenode > (length(genenodedepth)+1)/2) 
	{
	while(spnodedepth[father] <= depth)
	{
		if(father == root)
		{
			findnode <- father
			break
		}
		else
		{
			findnode <- father
			father <- sptree[father,1]
		}
	}
	}
	
	return (findnode)

}
			