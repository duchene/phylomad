`partition.tree` <-
function(tree,nspecies)
{
   partition<-matrix(0,(nspecies-2),nspecies)
   for(i in (nspecies+1):(2*nspecies-2)){
	group<-offspring.species(i,tree,nspecies)
	partition[(i-nspecies),group]<-1
   }
   partition<-cbind(partition,matrix(1,dim(partition)[1],1))
   return(partition)
}
