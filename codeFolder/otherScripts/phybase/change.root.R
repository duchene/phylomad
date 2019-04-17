`change.root` <-
function(nodematrix, newroot)
{

   root<-rootoftree(nodematrix)

 if(newroot != root){
    ##unroot the tree
    if(nodematrix[root,1] == -9)
	nodes<-unroottree(nodematrix)   else
	nodes<-nodematrix
    nspecies<-dim(nodes)[1]/2+1
    ##new root should be the internodes
    if(newroot<=nspecies)
	stop("new root should be the internodes!")

    anc<-ancestor(newroot,nodes)
    nanc<-length(anc)
    ##update old root
    a<-nodes[dim(nodes)[1],2:4]
    a[a==anc[nanc-1]]<-a[3]
    nodes[dim(nodes)[1],2:4]<-a

    ##reverse the ancestral history
    if(nanc==2){
	##reverse son, father, and branch length
	brlens<-nodes[anc[1],4]
	nodes[anc[1],4]<-anc[2]
	nodes[anc[2],1]<-anc[1]
	nodes[anc[2],4]<-brlens
    }
	
    if(nanc>2){
	##reverse son
	brlens<-nodes[anc[1],4]
	nodes[anc[1],4]<-anc[2]
	for(i in 2:(nanc-1)){
		if(nodes[anc[i],2]==anc[i-1])
			nodes[anc[i],2]<-anc[i+1]
		else
			nodes[anc[i],3]<-anc[i+1]
	}
	##reverse son and branch length
	nodes[anc[2:nanc],1]<-anc[1:(nanc-1)]
	nodes[anc[2:nanc],4]<-nodes[anc[1:(nanc-1)],4]
	nodes[anc[2],4]<-brlens
    }
    nodes[newroot,1]<--8
  }else{
    nodes<-nodematrix
  }
		
    z <- list(nodes = matrix(0, dim(nodes)[1], 5), rootnode=as.integer)

    z$nodes<-nodes
    z$rootnode<-newroot
    z	
}