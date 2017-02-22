get.biodivstat <- function(al){
	if(class(al) != "DNAbin"){ al <- as.character(as.DNAbin(al)) } else { al <- as.character(al) }
	div <- apply(al, 2, unique)
	div <- unlist(lapply(div, function(x) length(x[x != "-"])))
	meandiv <- mean(div[div != 0])
	return(meandiv)
}