get.model <- function(sdata, format = "phyllip"){
	  
	 candidates <- c("GTR+G", "GTR", "HKY+G", "HKY", "JC+G", "JC")

	 if(format == "phylip"){
                  data <- read.dna(sdata)
         } else if(format == "fasta"){
                  data <- read.dna(sdata, format = "fasta")
         } else if(format == "nexus"){
	   	  data <- as.DNAbin(read.nexus.data(sdata))
	 }

	 ICs <- modelTest(as.phyDat(data), model = c("GTR", "HKY", "JC"), G = T, I = F)
	 model <- candidates[which(ICs$BIC == min(ICs$BIC))]

	 print(ICs)

	 return(model)
	 
}