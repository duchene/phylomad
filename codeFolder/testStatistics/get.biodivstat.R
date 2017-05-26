get.biodivstat <- function(al){
	
	aadata <- class(al) == "AAbin"
        if(class(al) == "DNAbin" | class(al) == "AAbin") al <- as.character(al)

	div <- apply(al, 2, unique)
	div <- unlist(lapply(div, function(x) length(x[x != "-"])))
	meandiv <- mean(div[div != 0])
	return(meandiv)
}