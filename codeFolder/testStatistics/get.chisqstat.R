get.chisqstat <- function(al){

	if(class(al) != "DNAbin"){ al <- as.character(as.DNAbin(al)) } else { al <- as.character(al) }
        
	cont <- matrix(NA, nrow = nrow(al), ncol = 4)
	for(i in 1:nrow(al)) cont[i,] <- table(al[i,])[c("a", "c", "g", "t")]
	chisqstat <- chisq.test(cont)$statistic

        return(chisqstat)

}