`unroottree` <-
function(nodematrix)
{

    nodes<-nodematrix
    root<-rootoftree(nodematrix)
    newroot<-nodes[root,2]
    nodes[nodes[root,3],4]<-nodes[newroot,4]+nodes[nodes[root,3],4]
    nodes[nodes[root,3],1]<-newroot
    nodes[newroot,4]<-nodes[root,3]
    nodes[newroot,1]<--8
    unrootnodes<-nodes[1:(dim(nodematrix)[1]-1),]

    return(unrootnodes)
  
}

