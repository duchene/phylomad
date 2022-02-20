# The following samples a quartet of taxa given an edge. Assumes an unrooted tree.

sample.quartet <- function(tr, edge){
	node1 <- tr$edge[edge, 2]
	desc1 <- Descendants(tr, node = node1, "children")
	tax1 <- tr$tip.label[Descendants(tr, node = desc1[1], "tips")[[1]]]
	samp1 <- sample(tax1, 1)

	tax2 <- tr$tip.label[Descendants(tr, node = desc1[2], "tips")[[1]]]
	samp2 <- sample(tax2, 1)

	tax3 <- tr$tip.label[Descendants(tr, node = Siblings(tr, node1), "tips")[[1]]]
	samp3 <- sample(tax3, 1)
	tax4 <- tr$tip.label[which(!tr$tip.label %in% c(tax1, tax2, tax3))]
	samp4 <- sample(tax4, 1)
	quartsamp <- drop.tip(tr, tr$tip.label[which(!tr$tip.label %in% c(samp1, samp2, samp3, samp4))])
	return(quartsamp)
}