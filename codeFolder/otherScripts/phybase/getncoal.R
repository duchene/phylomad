'getncoal'<-
function(inode,sptree,nspecies, species.structure,coal,coalt)
{

	if(inode <= nspecies)
	{
		coalt[inode,1]<-sum(species.structure[inode,])
		coalt[inode,2]<-sum(coal[,1]==inode)
		coalt[inode,3]<-coalt[inode,1]-coalt[inode,2]
		return(coalt)
		
	}
	else
	{
		coalt<-getncoal(sptree[inode,2],sptree,nspecies,species.structure,coal,coalt)
		coalt<-getncoal(sptree[inode,3],sptree,nspecies,species.structure,coal,coalt)
		
		coalt[inode,1]<-coalt[sptree[inode,2],3]+coalt[sptree[inode,3],3]
		coalt[inode,2]<-sum(coal[,1]==inode)
		coalt[inode,3]<-coalt[inode,1]-coalt[inode,2]			
		return(coalt)
	}
}
