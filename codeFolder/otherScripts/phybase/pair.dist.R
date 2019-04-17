`pair.dist` <-
function(nodematrix,nspecies)
{
    dist<-matrix(0,nspecies,nspecies)
    for(i in 1:(nspecies-1)){
	for(j in (i+1):nspecies){
		dist[i,j]<-mrca.2nodes(i,j,nodematrix)$dist
        }
    }
    return(dist+t(dist))
}