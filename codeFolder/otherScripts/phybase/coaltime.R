`coaltime` <-
function(inode,jnode,nodematrix,nspecies)
{
    mrca<-mrca.2nodes(inode,jnode,nodematrix)$anc
    return(node.height(mrca,nodematrix,nspecies))

}

