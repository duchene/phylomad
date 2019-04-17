`treedist` <-
function(tree1,tree2)
{
    	nodes1<-tree1
    	nodes2<-tree2
    	if(!is.rootedtree(tree1))
    		nodes1<-root.tree(tree1,1)
    	if(!is.rootedtree(tree2))
    		nodes2<-root.tree(tree2,1)

    	nspecies<-(dim(nodes1)[1]+1)/2
    	partition1<-partition.tree(nodes1,nspecies)
    	partition2<-partition.tree(nodes2,nspecies)

    	dist<-dim(partition1)[1]

    	for(i in 1:dim(partition2)[1])
	{
		found<-FALSE
		for(j in 1:dim(partition1)[1])
		{
			if(sum(abs(partition2[i,]-partition1[j,]))==0)
			{
				dist<-dist-1
				found<-TRUE
				break
			}
    		}
		if(!found)
		dist<-dist+1
    	}
    	return(dist)
}