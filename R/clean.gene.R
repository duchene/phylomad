clean.gene <- function(x){
	   gene <- as.character(x)
	   goodcols <- vector()
	   for(i in 1:ncol(gene)){
	   	 if(length(grep("a|c|g|t", gene[, i])) == nrow(gene)) goodcols <- c(goodcols, i)
	   }
	   if(length(goodcols) == 0 || length(goodcols) == 1) stop("Too many sites have missing data")
	   gene <- gene[, goodcols]
	   print(length(goodcols))
	   return(as.DNAbin(gene))
}