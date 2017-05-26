get.chisqstat <- function(al){
	
	aadata <- class(al) == "AAbin"
	if(class(al) == "DNAbin" | class(al) == "AAbin") al <- as.character(al)
        
	al <- tolower(al)
	
	if(aadata){
		cont <- matrix(NA, nrow = nrow(al), ncol = 20)
		for(i in 1:nrow(al)) cont[i,] <- table(al[i,])[c("c", "d", "s", "q", "k", "i", "p", "t", "f", "n", "g", "h", "l", "r", "w", "a", "v", "e", "y", "m")]
	} else {
		cont <- matrix(NA, nrow = nrow(al), ncol = 4)
	        for(i in 1:nrow(al)) cont[i,] <- table(al[i,])[c("a", "c", "g", "t")]
	}
	
	chisqstat <- chisq.test(cont)$statistic

        return(chisqstat)

}