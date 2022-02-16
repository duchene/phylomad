clean.gene <- function(sdata, format = "phylip", aadata = F, clean = T, clean.by.edge = T){
	   if(format == "phylip"){
                  if(aadata) data <- as.AAbin(read.aa(sdata)) else data <- read.dna(sdata)
           } else if(format == "fasta"){
                  if(aadata) data <- read.aa(sdata, format = "fasta") else data <- read.dna(sdata, format = "fasta")
           } else if(format == "nexus"){
		  if(aadata) data <- as.AAbin(read.nexus.data(sdata)) else data <- as.matrix(as.DNAbin(read.nexus.data(sdata)))
	   }
	   
	   
	   gene <- as.character(as.matrix(data))
	   
	   goodtax <- vector()
	   goodcol <- vector()
	   if(clean.by.edge){
		if(aadata){
			for(i in 1:nrow(gene)) if(length(which(gene[i,] %in% c("?", "-", "O", "o", "X", "x"))) != ncol(gene)) goodtax <- c(goodtax, i)
			for(i in 1:ncol(gene)) if(length(which(gene[,i] %in% c("?", "-", "O", "o", "X", "x"))) != nrow(gene) & length(grep("c|d|s|q|k|i|p|t|f|n|g|h|l|r|w|a|v|e|y|m|C|D|S|Q|K|I|P|T|F|N|G|H|L|R|W|A|V|E|Y|M", gene[, i])) >= (nrow(gene) * 0.05)) goodcol <- c(goodcol, i)
	   	} else {
			for(i in 1:nrow(gene)) if(length(which(gene[i,] %in% c("N", "n", "?", "-", "O", "o", "X", "x"))) != ncol(gene)) goodtax <- c(goodtax, i)
			for(i in 1:ncol(gene)) if(length(which(gene[,i] %in% c("N", "n", "?", "-", "O", "o", "X", "x"))) != nrow(gene) & length(grep("a|c|g|t|A|C|G|T", gene[, i])) >= (nrow(gene) * 0.05)) goodcol <- c(goodcol, i)
		}
		if(length(goodtax) < 3) stop("Too many taxa have fully missing data")
		if(length(goodcol) < 10) stop("Too many columns have fully missing data")
	   	gene <- gene[goodtax, goodcol]
	   }
	   
	   if(!clean & aadata) return(as.AAbin(gene)) else if(!clean) return(as.DNAbin(gene))
	   
	   goodcols <- vector()
	   if(aadata){
		for(i in 1:ncol(gene)) if(length(grep("c|d|s|q|k|i|p|t|f|n|g|h|l|r|w|a|v|e|y|m|C|D|S|Q|K|I|P|T|F|N|G|H|L|R|W|A|V|E|Y|M", gene[, i])) >= (nrow(gene) * 0.95)) goodcols <- c(goodcols, i)
	   } else {
		for(i in 1:ncol(gene)) if(length(grep("a|c|g|t|A|C|G|T", gene[, i])) >= (nrow(gene) * 0.95)) goodcols <- c(goodcols, i)
	   }
	   if(length(goodcols) == 0 || length(goodcols) == 1) stop("Too many sites have missing data")
	   if(length(goodcols) < (ncol(gene) * 0.1)) stop("Less than 10% of sites can be used for testing saturation")
	   gene <- gene[, goodcols]
	   print(paste("Gene was left with", length(goodcols), "sites."))
	   
	   if(aadata) return(as.AAbin(gene)) else return(as.DNAbin(gene))
}