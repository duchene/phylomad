# Added by Brian O'Meara
phybase2phylo <- function(x) {
	if(grepl('#', x)) {
		x<-gsub("#\\.*\\d*\\.*\\d*","",x)
	}
	return(read.tree(text=x))
}
