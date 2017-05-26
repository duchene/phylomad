clean.gene <- function(sdata, format = "phylip", aadata = F){
	   if(format == "phylip"){
                  if(aadata) data <- as.AAbin(read.aa(sdata)) else data <- read.dna(sdata)
           } else if(format == "fasta"){
                  if(aadata) data <- read.aa(sdata, format = "fasta") else data <- read.dna(sdata, format = "fasta")
           } else if(format == "nexus"){
		  if(aadata) data <- as.AAbin(read.nexus.data(sdata)) else data <- as.DNAbin(read.nexus.data(sdata))
	   }
	   
	   gene <- as.character(as.matrix(data))
	   goodcols <- vector()
	   if(aadata){
		for(i in 1:ncol(gene)) if(length(grep("c|d|s|q|k|i|p|t|f|n|g|h|l|r|w|a|v|e|y|m|C|D|S|Q|K|I|P|T|F|N|G|H|L|R|W|A|V|E|Y|M", gene[, i])) == nrow(gene)) goodcols <- c(goodcols, i)
	   } else {
		for(i in 1:ncol(gene)) if(length(grep("a|c|g|t|A|C|G|T", gene[, i])) == nrow(gene)) goodcols <- c(goodcols, i)
	   }
	   if(length(goodcols) == 0 || length(goodcols) == 1) stop("Too many sites have missing data")
	   gene <- gene[, goodcols]
	   print(paste("Gene was left with", length(goodcols), "sites."))
	   if(aadata) return(as.AAbin(gene)) else return(as.DNAbin(gene))
}