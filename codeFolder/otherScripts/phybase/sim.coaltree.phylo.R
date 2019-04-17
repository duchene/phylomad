sim.coaltree.phylo <- function(phy, pop.size=1, nsamples=1) {
	phybase.tree <- write.tree(phy)
	spname <- species.name(phybase.tree)
	nodematrix <- read.tree.nodes(phybase.tree, spname)$nodes
	nodematrix[,5] <- pop.size
	root.node <- rootoftree(nodematrix)
	sequences <- rep(nsamples, length(spname))
	return(phybase2phylo(sim.coaltree.sp(root.node, nodematrix, nspecies=length(spname), sequences, spname)$gt))
}
