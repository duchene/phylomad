clean.gene <- function(sdata, format = "phylip"){
	   if(format == "phylip"){
                  data <- read.dna(sdata)
           } else if(format == "fasta"){
                  data <- read.dna(sdata, format = "fasta")
           } else if(format == "nexus"){
		  data <- as.DNAbin(read.nexus.data(sdata))
	   }
	   gene <- as.character(as.matrix(data))
	   goodcols <- vector()
	   for(i in 1:ncol(gene)){
	   	 if(length(grep("a|c|g|t", gene[, i])) == nrow(gene)) goodcols <- c(goodcols, i)
	   }
	   if(length(goodcols) == 0 || length(goodcols) == 1) stop("Too many sites have missing data")
	   gene <- gene[, goodcols]
	   print(paste("Gene was left with", length(goodcols), "sites."))
	   return(as.DNAbin(gene))
}