`sim.dna` <-
function(nodematrix,seqlength,model,kappa=2,rate=c(1,1,1,1,1,1),frequency=c(1/4,1/4,1/4,1/4))
{

	############################################################################################
	# 	model = 1 : JC model
	# 	model = 2 :  H2P model, alpha and beta (two parameters).
	# 	model = 3 : HKY model, alpha, beta, frequencies (five parameters)
	# 	model = 4 : GTR model, six rates, frequencies (nine parameters)
	#
	#	rate[1] = A <-> C rate
	#     rate[2] = A <-> G rate
	#     rate[3] = A <-> T rate
	#     rate[4] = C <-> G rate
	#     rate[5] = C <-> T rate
	#     rate[6] = G <-> T rate
	#
	#     code for constructing rate matrix is copied from Mrbayes (Huelsenbeck and Ronquist) 
	############################################################################################

 	
	if(model==1)
	{
		rateValues<-rep(1,6)
		bs<-rep(1/4,4)
	}
	else if(model == 2)
	{
		rateValues<-rep(1,6)
		rateValues[c(2,5)]<-kappa
		bs<-rep(1/4,4)
	}
	else if(model == 3)
	{
		rateValues<-rep(1,6)
		rateValues[c(2,5)]<-kappa
		bs<-frequency
	}
	else 
	{
		rateValues <- rate
		bs <- frequency
	}

	a<-matrix(0,4,4)
	scaler <- 0.0
	for (i in 1:3)
	{
	 	for (j in (i+1):4)
		{
				{if (i == 1 && j == 2)
					mult <- rateValues[1]
				else if (i == 1 && j == 3)
					mult <- rateValues[2]
				else if (i == 1 && j == 4)
					mult <- rateValues[3]
				else if (i == 2 && j == 3)
					mult <- rateValues[4]
				else if (i == 2 && j == 4)
					mult <- rateValues[5]
				else if (i == 3 && j == 4)
					mult <- rateValues[6]}
				a[i,j] <- bs[j] * mult
				a[j,i] <- bs[i] * mult
				a[i,i] <- a[i,i] - a[i,j]
				a[j,j] <- a[j,j] - a[j,i]
				scaler <- scaler + bs[i] * a[i,j]
				scaler <- scaler + bs[j] * a[j,i]
		}
	}
			
	#rescale Q matrix 
	scaler = 1.0 / scaler;
	for (i in 1:4)
 		for (j in 1:4)
			a[i,j] = a[i,j] * scaler
		

    	root<-rootoftree(nodematrix)
    	{if(nodematrix[root,1] == -8){
    		sonnodes<-nodematrix[root,2:4]
		nsonnodes<-3
  		nsequence<-(dim(nodematrix)[1]+2)/2
    	} 
	else{
		sonnodes<-nodematrix[root,2:3]
		nsonnodes<-2
		nsequence<-(dim(nodematrix)[1]+1)/2
    	}}
    ##nucleotides at root node follow the equilibrium distribution defined by parameter frequency.
    	dna<-matrix(0,nrow=dim(nodematrix)[1],ncol=seqlength) 
    	dna[root, ] <- sample(1:4,size=seqlength,replace=TRUE,prob=frequency)
    	tranprob<-matrix(0,seqlength,4)
	

    	##nucleotides at other nodes 
	while(nsonnodes>0){
		for(i in 1:nsonnodes){
			father<-nodematrix[sonnodes[i],1]
			dna[sonnodes[i],]<-simnucleotide(dna[father,],nodematrix[sonnodes[i],4], a)
     		}
		sonnodes<-nodematrix[sonnodes,2:3]
		sonnodes<-sonnodes[sonnodes>0]
		{if(length(sonnodes)>0)
			nsonnodes<-length(sonnodes)
		 else
			nsonnodes<-0
		}
    	}
    	return(dna[1:nsequence,])	
}