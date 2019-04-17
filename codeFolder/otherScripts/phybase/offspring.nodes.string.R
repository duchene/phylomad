`offspring.nodes.string` <-
function(inode,nodematrix,nspecies)
{
    if(inode<1 | inode>dim(nodematrix)[1]){
a<-paste("The node number should be between 1 and",dim(nodematrix)[1])
stop(a)
    } 
    if(inode<=nspecies)
return(paste(inode))
    if(inode>nspecies){
son1<-nodematrix[inode,2]
  son2<-nodematrix[inode,3]
str1<-offspring.nodes.string(son1,nodematrix,nspecies)
        str2<-offspring.nodes.string(son2,nodematrix,nspecies)
str<-paste(inode)
str<-paste(str,str1)
        return(paste(str,str2))
    }
}

