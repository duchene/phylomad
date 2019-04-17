`popsize` <-
function(inode,jnode,nodematrix)
{
    mrca<-mrca.2nodes(inode,jnode,nodematrix)$anc
    return(nodematrix[mrca,5])

}

