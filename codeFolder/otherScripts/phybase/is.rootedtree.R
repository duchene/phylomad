`is.rootedtree` <-
function(tree)
{
	if(is.character(tree))
	{
		string <- unlist(strsplit(tree, NULL))
    		leftpar<-which(string=="(")  
    		rightpar<-which(string==")")
		comma<-which(string==",") 
    		if(length(leftpar) != length(leftpar))
     			stop("The number of left parenthesis is NOT equal to the number of right parenthesis")
    		if(length(leftpar) == length(comma))
			rooted<-TRUE
    		else if(length(leftpar) == (length(comma)-1))
			rooted<-FALSE
    		else
			stop("The number of comma in the tree string is wrong!")
    	} 
    	else
	{
		rooted<-TRUE
		if(tree[rootoftree(tree),1] == -8)
			rooted<-FALSE
    	}
    	return(rooted)
}

