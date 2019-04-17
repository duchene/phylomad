steac.sptree<-function(trees,speciesname,taxaname,species.structure,outgroup,method="nj")
{
    	ntree<-length(trees)
    	nspecies<-length(speciesname)
    	ntax<-length(taxaname)
    	dist <- matrix(0, nrow=ntree, ncol=ntax*ntax)


    	for(i in 1:ntree)
	{
  		ranktree1 <- read.tree.nodes(trees[i])
		thistreetaxa<-ranktree1$names
		ntaxaofthistree<-length(thistreetaxa)
		ranktree<-ranktree1$nodes
		thistreenode<-rep(-1,ntaxaofthistree)

		#find the missing taxa
		for(j in 1:ntaxaofthistree)
		{
			thistreenode[j]<-which(taxaname == thistreetaxa[j])
			if(length(thistreenode[j])==0)
			{
			print(paste("wrong taxaname",thistreetaxa[j],"in gene",i))
			return(0)
			}
		}


		#if(!is.rootedtree(ranktree))
		#{
		#	root<-which(thistreetaxa==outgroup)
		#	if(length(root) == 0)
		#		warnings(paste("outgroup is missing at tree",i))
		#	else
		#		ranktree<-root.tree(ranktree,root)
		#}
	
		#estimate the coalescence times
  		#treestr<-write.subtree(rootoftree(ranktree),ranktree,ntaxaofthistree,rootoftree(ranktree))
		#tree<-node2name(treestr, thistreetaxa)
		#noclocktree<-read.tree(text=tree)
		#clocktree<-write.tree(chronogram(noclocktree))
		#clocktree<-read.tree.nodes(clocktree,thistreetaxa)$nodes
		#coaltimes<-rep(0,ntaxaofthistree-1)
		#for(j in (ntaxaofthistree+1):(2*ntaxaofthistree-1))
		#	coaltimes[j-ntaxaofthistree]<-node.height(j,clocktree,ntaxaofthistree)

		dist1 <- matrix(0, ntax, ntax)
		for (k in 1:(ntaxaofthistree - 1)) {
        		for (j in (k + 1):ntaxaofthistree) {
            		dist1[thistreenode[k], thistreenode[j]] <- mrca.2nodes(k, j, ranktree)$dist
        		}
    		} 
   
		dist[i,]<-as.numeric(dist1+t(dist1))
    	}
    	dist[dist==0]<-NA
    	dist2<-matrix(apply(dist,2,mean,na.rm=TRUE),ntax,ntax)
    	diag(dist2)<-0
    	if(sum(is.nan(dist2))>0)
    	{
		print("missing species!")
		dist2[is.nan(dist2)]<-10000
    	}
    	dist3 <- pair.dist.mulseq(dist2, species.structure)
    	dist <- dist3
    	{
		if(method == "upgma") sptree2 <- upgma(dist,speciesname)$treestr
     		if(method == "nj") 
		{
		sptree<-write.tree(nj(dist))
		sptree1<-read.tree.nodes(sptree,speciesname)$nodes
		root<-which(taxaname==outgroup)
		root<-species.structure[,root]
		root<-which(root==1)
		sptree2<-root.tree(sptree1,root)
		a<-sptree2[,4]
		b<-which(a<0)
		if(length(b) > 1)
			sptree2[b[1:(length(b)-1)],4]<-0
		sptree2<-write.subtree(rootoftree(sptree2),sptree2,speciesname,rootoftree(sptree2))
		}
		
	     	#else print("choose either upgma or nj")
    	}

    	sptree<-sptree2
    	return(sptree)
}
