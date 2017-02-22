get.model <- function(sdata, format = "phyllip"){
	  
	 candidates <- c("GTR", "GTR+G")

	 if(format == "phyllip"){
                  data <- read.dna(sdata)
         } else if(format == "fasta"){
                  data <- read.dna(sdata, format = "fasta")
         }	  

	 ICs <- modelTest(as.phyDat(data), model = "GTR", I = F)
	 model <- candidates[which(ICs$BIC == min(ICs$BIC))]

	 return(model)
	 
}