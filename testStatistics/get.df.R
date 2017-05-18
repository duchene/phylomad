get.df <- 
function(tre){
	brlens <- as.numeric(scale(tre$edge.length, scale = F, center = T))
	df <- sum(brlens[which(tre$edge[,2] <= Ntip(tre))]) - sum(brlens)
	return(df)
}

