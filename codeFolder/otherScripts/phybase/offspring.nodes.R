`offspring.nodes` <-
function(inode,nodematrix,nspecies)
{
    	string<-offspring.nodes.string(inode,nodematrix,nspecies)
    	offspringnodes<-as.numeric(unlist(strsplit(string,split=" ")))
    	if(nodematrix[inode,1] == -8)
		offspringnodes<-dim(nodematrix)[1]:1
    	return(as.numeric(unlist(strsplit(string,split=" "))))
}

