`mrca.nodes` <-
function(nodevector,nodematrix)
{
	nnodes<-length(nodevector)
    	if(nnodes<2)
		stop("The number of nodes should be at least 2!")

    	mrca<-mrca.2nodes(nodevector[1],nodevector[2],nodematrix)$anc
    	for(i in 3:nnodes)
		mrca<-mrca.2nodes(mrca,nodevector[i],nodematrix)$anc

    	return(mrca)
}

