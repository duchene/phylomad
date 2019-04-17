`postdist.tree` <-
function(trees,name="")
{
    ntree<-length(trees)
    if(length(name)>0)
    	name<-sort(species.name(trees[1]))
    nspecies<-length(name)
    treestr<-rep("",ntree)
    weight<-rep(0,ntree)
    tree<-trees
    index<-1

    if(length(grep("[a-z]",name[1],ignore.case=TRUE))>0){
    while(length(tree)>0){
  	ntree<-length(tree)
	ignoretrees<-rep(0,ntree)
	for(i in 1:length(tree)){
		treestr[index]<-tree[1]
		nodematrix1<-read.tree.nodes(tree[1],name)$nodes
		nodematrix2<-read.tree.nodes(tree[i],name)$nodes
		dist<-treedist(nodematrix1,nodematrix2)
		if(dist == 0){
			weight[index]<-weight[index]+1
			ignoretrees[i]<-1
  		}
	}
	index<-index+1
 	notignoretrees<-(1-ignoretrees)*(1:ntree)
	tree<-tree[notignoretrees]
    }
    } else{
    while(length(tree)>0){
  	ntree<-length(tree)
	ignoretrees<-rep(0,ntree)
	for(i in 1:length(tree)){
		treestr[index]<-tree[1]
		nodematrix1<-read.tree.nodes(tree[1])$nodes
		nodematrix2<-read.tree.nodes(tree[i])$nodes
		dist<-treedist(nodematrix1,nodematrix2)
		if(dist == 0){
			weight[index]<-weight[index]+1
			ignoretrees[i]<-1
  		}
	}
	index<-index+1
 	notignoretrees<-(1-ignoretrees)*(1:ntree)
	tree<-tree[notignoretrees]
    }
    }

    z<-list(tree=as.vector,prob=as.vector)
    z$tree<-treestr[treestr!=""]
    z$prob<-weight[weight>0]/length(trees)
    z


}
