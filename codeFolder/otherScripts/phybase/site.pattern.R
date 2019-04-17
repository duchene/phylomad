'site.pattern' <-
function(seq)
{
	numsite<-dim(seq)[2]
	numtaxa<-dim(seq)[1]
	
	result<-matrix(0,nrow=numtaxa,ncol=numsite)

	for(i in 2:numtaxa)
	{
  		b<-rep(TRUE,numsite)
  		for(j in 1:(i-1))
   		{
			b[seq[i,]==seq[j,]]<-FALSE			
			result[i,b]<-result[j,b]+1
   		}
	}

	result<-sortmat(t(result),1:numtaxa)
	index<-1
    	number<-1
	numpattern<-rep(-1,numsite)
	pattern<-matrix(-1,nrow=numsite,ncol=numtaxa)

	for(i in 2:numsite)
	{
		
		if(sum(abs(result[i,]-result[i-1,])) == 0)
		{
			number<-number + 1
		}
		else
		{
			numpattern[index]<-number
			pattern[index,]<-result[i-1,]
			index<-index+1
			number<-1
		}
	}
	numpattern[index]<-number
	pattern[index,]<-result[numsite,]
	
    	return(cbind(pattern[1:index,],numpattern[1:index])) 
}