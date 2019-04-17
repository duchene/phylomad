'loglikeSP'<-
function(gtree,sptree,taxaname, spname, species.structure, strict=T)
{
	ntax<-dim(species.structure)[2]
	nspecies<-dim(species.structure)[1]
	ntree<-length(gtree)
	if(!is.matrix(sptree))sptree<-read.tree.nodes(sptree,spname)$nodes
	rootofsptree<-rootoftree(sptree)
	ab<-matrix(0,nrow=2*nspecies-1,ncol=2)
		
	for(k in 1:ntree)
	{
		genetree<-read.tree.nodes(gtree[k],taxaname)$nodes
		b<-rep(-1,2*nspecies-1)

		coal<-getcoaltime(genetree,sptree,ntax,nspecies,species.structure)
		coalt<-matrix(-1,nrow=2*nspecies-1,ncol=3)
		coalt<-getncoal(rootofsptree,sptree,nspecies,species.structure,coal,coalt)

		node<-which(coalt[,1]>1)
		for(i in 1:length(node))
		{
			height<-sptree[node[i],4]
			time<-rep(height,coalt[node[i],2]+1)
			if(coalt[node[i],2]>0) time[1:coalt[node[i],2]]<-sort(coal[which(coal==node[i]),2])
			time1<-rep(0,length(time))
			if(length(time)>1) time1[2:length(time)]<-time[1:(length(time)-1)]
			time<-time-time1
			n<-coalt[node[i],1]:(coalt[node[i],3])
			b[node[i]]<-sum(n*(n-1)*time)
		}

		ab[,1]<-ab[,1]+coalt[,2]
		ab[,2]<-ab[,2]+b
	}

	ab<-cbind(ab,sptree[,5])
	ab<-ab[which(ab[,2]>-1),]
	if(strict){if(sum(ab[,1]) != (ntax-1)*ntree) stop("the total number of coalescence is not equal to ntax-1")}
	loglike<-sum(ab[,1]*log(2/ab[,3])-(ab[,2]/ab[,3]))

	return(loglike)
}	
