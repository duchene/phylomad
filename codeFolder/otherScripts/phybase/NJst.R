'NJst'<-
function(genetrees, taxaname, spname, species.structure)
{

	ntree<-length(genetrees)
	ntaxa<-length(taxaname)
	dist <- matrix(0, nrow = ntree, ncol = ntaxa * ntaxa)
	
	for(i in 1:ntree)
	{
		genetree1 <- read.tree.nodes(genetrees[i])
        	thistreetaxa <- genetree1$names
        	ntaxaofthistree <- length(thistreetaxa)
        	thistreenode <- rep(-1, ntaxaofthistree)
		dist1<-matrix(0,ntaxa,ntaxa)
        	for (j in 1:ntaxaofthistree)
		{
            		thistreenode[j] <- which(taxaname == thistreetaxa[j])
            		if (length(thistreenode[j]) == 0)
			{
                		print(paste("wrong taxaname", thistreetaxa[j],"in gene", i))
                		return(0)
            		}
        	}
		dist1[thistreenode, thistreenode]<-nancdist(genetrees[i],thistreetaxa)$dist
		dist[i,]<-as.numeric(dist1)
	}

	dist[dist == 0] <- NA
    	dist2 <- matrix(apply(dist, 2, mean, na.rm = TRUE), ntaxa, ntaxa)
    	diag(dist2) <- 0
    	if (sum(is.nan(dist2)) > 0)
	{
        	print("missing species!")
        	dist2[is.nan(dist2)] <- 10000
    	}
    	speciesdistance <- pair.dist.mulseq(dist2, species.structure)

	tree<-write.tree(nj(speciesdistance))
	node2name(tree,name=spname)
}
