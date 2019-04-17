triplenumber <- function(dna)
{
	ntaxa<-dim(dna)[1]
	ntriple<-ntaxa*(ntaxa-1)*(ntaxa-2)/6
	triple<-matrix(0,nrow=ntriple,ncol=5)

	for(i in 1:(ntaxa-2))
		for(j in (i+1):(ntaxa-1))
			for(k in (j+1):ntaxa)
			{
				seq<-dna[c(i,j,k),]
				triple[i+j+k-5,]<-site.pattern(seq)[,4]
  			}

	triple
			
}