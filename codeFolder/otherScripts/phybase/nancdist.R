'nancdist'<-
function(tree, taxaname)
{
	ntaxa<-length(taxaname)
	nodematrix<-read.tree.nodes(tree,taxaname)$nodes
	if(is.rootedtree(nodematrix)) nodematrix<-unroottree(nodematrix)
	dist<-matrix(0, ntaxa,ntaxa)
	for(i in 1:(ntaxa-1))
		for(j in (i+1):ntaxa)
		{
		anc1<-ancestor(i,nodematrix)
		anc2<-ancestor(j,nodematrix)
		n<-sum(which(t(matrix(rep(anc1,length(anc2)),ncol=length(anc2)))-anc2==0, arr.ind=TRUE)[1,])-3
		if(n==-1) n<-0
		dist[i,j]<-n
		}
	dist<-dist+t(dist)
	z<-list(dist=as.matrix, taxaname=as.vector)
	z$dist<-dist
	z$taxaname<-taxaname
	z
}

