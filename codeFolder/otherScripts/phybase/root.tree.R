`root.tree` <-
function(nodematrix,outgroup)
{


    node<-change.root(nodematrix,nodematrix[outgroup,1])
    nodes<-node$nodes
    oldroot<-node$rootnode

    rootnode<-matrix(-9,1,dim(nodematrix)[2])
    rootnode[1,2]<-outgroup
    rootnode[1,3]<-oldroot
    

    if(nodes[oldroot,2]==outgroup)
	nodes[oldroot,2]<-nodes[oldroot,4] 
    if(nodes[oldroot,3]==outgroup)
	nodes[oldroot,3]<-nodes[oldroot,4] 

    nodes[oldroot,4]<-nodes[outgroup,4]/2
    nodes[outgroup,4]<-nodes[outgroup,4]/2

    nodes[oldroot,1]<-dim(nodes)[1]+1
    nodes[outgroup,1]<-dim(nodes)[1]+1

    roottree<-rbind(nodes,rootnode)
   
    return(roottree)

}

