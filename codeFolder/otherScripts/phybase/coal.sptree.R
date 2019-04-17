`coal.sptree` <-
function(trees,speciesname,nspecies,outgroup=1)
{

    ntree<-length(trees)
    dist<-matrix(0,nspecies,nspecies)
    totallength<-0

    for(i in 1:ntree){
	treenode<-read.tree.nodes(trees[i],speciesname)$nodes
	if(!is.rootedtree(treenode))
		treenode<-root.tree(treenode,outgroup)
	clocktree<-noclock2clock(2*nspecies-1, treenode, nspecies)
	dist1<-pair.dist(clocktree,nspecies)
	dist<-dist+dist1/sum(dist1)
 	totallength<-totallength+sum(dist1)
    }
    dist<-dist*totallength/(ntree*ntree)
    sptree<-upgma(dist,speciesname)
    z<-list(nodes=as.matrix,treestr="",name="")
    z$nodes<-sptree$nodes
    z$treestr<-sptree$treestr
    z$name<-sptree$speciesname
    z
}