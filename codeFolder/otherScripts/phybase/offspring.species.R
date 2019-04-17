`offspring.species` <-
function(inode,nodematrix,nspecies)
{   
    offspring<-offspring.nodes(inode,nodematrix,nspecies)
    return(offspring[offspring<=nspecies])
}

