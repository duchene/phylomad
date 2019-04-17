"populationMutation"<-
function (sptree, spnodedepth, genetree, genenodedepth, speciesmatrix)
{  
	index<-1

	nspecies<-(dim(sptree)[1]+1)/2
	ntaxa<-(dim(genetree)[1]+1)/2
	genetreenodes <- rep(-1,2*ntaxa)

	for(i in 1:nspecies)
      	{
		seq<-(1:ntaxa)*(speciesmatrix[i,])
		seq<-seq[seq>0]
        	for(inode in 1:length(seq))
        	{        
        		inodegene <- seq[inode];
        		stop=0;
        		while(inodegene != dim(genetree)[1])
			{
				#check if the node is already taken care of
				for(k in 1:index)
              				if(inodegene == genetreenodes[k]) 
					{
						stop<-1
						break
					}
				if(stop == 1) 	break

				#change the branch length of node p				
				genetree[inodegene,4] <- ChangeBrlen(sptree, spnodedepth, i, genetree, genenodedepth, inodegene)

				#copy p to genetreenode
				genetreenodes[index] <- inodegene
                 		index<-index+1

				#reset inodegene
				inodegene<-genetree[inodegene,1]
			}
		}
	}
	genetree[,4]<-round(genetree[,4],6)
	return (genetree)

}