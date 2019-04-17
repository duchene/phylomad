`pair.dist.mulseq` <-
function(dist,species.structure)
{
    for(i in 1:dim(species.structure)[1]){
	species.structure[i,]<-species.structure[i,]/sum(species.structure[i,])
    }

    dis<-round((species.structure)%*%dist%*%t(species.structure),8)
    diag(dis)<-0
    dis
}

