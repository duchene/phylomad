get.unconstrained <- function(al){
	
	aadata <- class(al) == "AAbin"
        if(class(al) == "DNAbin" | class(al) == "AAbin") al <- as.character(al)	  
	
	alstrs <- apply(al, 2, paste0, collapse = "")
	strstab <- as.numeric(table(alstrs))
	nsites <- ncol(al)
	uncL <- sum(sapply(strstab, function(x) x*log(x)))
	uncL <- uncL - (nsites*log(nsites))
	return(uncL)
}