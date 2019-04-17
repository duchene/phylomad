'write.dna'<-
function (sequence, name, file = "", format="nexus", program="mrbayes",partition=matrix(0,ncol=2,nrow=1), clock=0, popmupr=0, ngen=1000000,nrun=1,nchain=1,samplefreq=100,taxa=as.vector,burnin=1000,gamma="(3,0.02)", outgroup=1,outfile="",append = FALSE) 
{
	#the sequence matrix: row represents taxa, column represents sites, very important!
    	ntax <- dim(sequence)[1]
    	nchar <- dim(sequence)[2]
	options(scipen=3)

	if(nchar < ntax)
		warnings("check the sequence matrix to make sure that row is taxa and column is sites")

    	add <- nchar%%100
	if(add != 0)
		add <- 100-add

      totalsequence<-cbind(sequence,matrix("",ncol=add,nrow=ntax))
	newseq<-matrix("",nrow=(nchar+add)*ntax/100,ncol=2)

	for(i in 1:((nchar+add)/100))
	{
		newseq[((i-1)*ntax+1):(i*ntax),]<-cbind(name,apply(totalsequence[,((i-1)*100+1):(i*100)],1,paste,collapse="",sep=""))
	}
	newseq<-cbind(newseq,rep("\n",dim(newseq)[1]))
	newseq[,1]<-paste(newseq[,1],"        ")
	
	{if(tolower(format) == "nexus")
	{
    		cat("#NEXUS\nBEGIN DATA;\nDIMENSIONS  NTAX=",ntax, "NCHAR=",nchar,";\nFORMAT DATATYPE=DNA  INTERLEAVE=YES MISSING=? GAP=- ;\nMATRIX", "\n",t(newseq),";\nEND;\n\n\n",file=file,append=append)
	}
	else
	{
		if(dim(newseq)[1] > ntax)
			newseq[(ntax+1):dim(newseq)[1],1] <- ""
		a<-rep("",dim(newseq)[1])
		b<-(1:(dim(newseq)[1]/ntax))*ntax
		a[b]<-"\n"
		string<-paste(newseq[,1],newseq[,2], a, sep="")
		title<-paste(ntax,nchar)
		
		write.table(title,file=file,append=append,row.names=FALSE, col.names=FALSE, quote=FALSE)
		write.table(string,file=file,append=TRUE,row.names=FALSE, col.names=FALSE, quote=FALSE)

		return (1)
	}}
	
    	#print out commands for variaty of methods
    	if(tolower(program) == "best")
    	{
		cat("BEGIN mrbayes;\n\tset autoclose=yes nowarn=yes;\n\toutgroup",outgroup,";\n",file=file,append=TRUE)
     		if(dim(partition)[1]>1)
    		{
			gene <- paste("gene",1:dim(partition)[1],sep="")
			string <- paste("\tcharset",gene,"=",partition[,1],"-",partition[,2],";\n")
    			cat(string, "\tpartition currentpartition = ", dim(partition)[1],":", file=file,append=TRUE)

    			string <- paste(gene,",") 
			string[length(string)] <- paste(gene[length(gene)],";\n")
    			cat(string,"\tset partition = currentPartition;\n",file=file,append=TRUE)
    		}

    		if(length(taxa)>1)
    		{
			nspecies<-length(taxa)
			spname<-paste("S",1:nspecies,sep="")
			speciesmatrix<-matrix(0,nrow=nspecies,ncol=ntax)
			index<-1	
			for(i in 1:length(taxa))
			{
				for(j in 1:taxa[i])
				{
					speciesmatrix[i,index]<-1
					index<-index+1
				}
			}
			for(i in 1:nspecies)
			{
				a<-(1:ntax)*speciesmatrix[i,]
				a<-a[a>0]
				cat("\ttaxset",spname[i],"=",a,";\n",file=file,append=TRUE)
			}
    		}
		string <- unlist(strsplit(file, NULL))
		a<-which(string=="/")
		if(length(a) == 0)
			spfile <- file
		else
			spfile<-paste(string[(a[length(a)]+1):length(string)],collapse="")
		spfile<-paste(spfile,".sptree",sep="")
		
		if(clock == 1)
			cat("\tprset   brlenspr=clock:coalescence"," popmupr=", popmupr," thetapr=invgamma", gamma," BEST=1;\n\tunlink  topology=(all) brlens=(all) ;\n\tmcmc ngen=",ngen," nrun = ",nrun," nchain = ",nchain," samplefreq=",samplefreq," ;\n\tsumt nrun=", nrun, "burnin = ",burnin, " filename=",spfile,"contype=allcompat;\n\tquit;\nEND;",file=file,append=TRUE)
           else
            	cat("\tprset "," popmupr=", popmupr," thetapr=invgamma", gamma," BEST=1;\n\tunlink  topology=(all) brlens=(all) ;\n\tmcmc ngen=",ngen," nrun = ",nrun," nchain = ",nchain," samplefreq=",samplefreq," ;\n\tsumt nrun=", nrun, "burnin = ",burnin, " filename=",spfile,"contype=allcompat;\n\tquit;\nEND;",file=file,append=TRUE)

    	}
    	if(program=="mrbayes")
    	{
		string <- unlist(strsplit(file, NULL))
		a<-which(string=="/")
		if(length(a) == 0)
			spfile <- file
		else
			spfile<-paste(string[(a[length(a)]+1):length(string)],collapse="")
		cat("begin mrbayes;\n\tset autoclose=yes nowarn=yes;\n\toutgroup",outgroup,";\n\tmcmc ngen=",ngen," nrun = ",nrun," nchain = ",nchain," samplefreq=",samplefreq," ;\n\tsumt nrun=", nrun, "burnin = ",burnin, " filename=",spfile,"contype=allcompat;\n\tquit;\nEND;",file=file,append=TRUE)
    	}
    	if(program=="paup")
    	{
		cat("begin paup;\nset defaultmode=yes;\noutgroup", outgroup,";\nlscore 1/nst=1;\nnj;\nroottrees root=outgroup;\nsavetrees file=",outfile," brlens=yes append=yes;\nset criterion=likelihood;\nlset nst=1 clock=no;\nalltrees;\nroottrees root=outgroup;\nsavetrees file=",outfile," brlens=yes append=yes;\nquit;\nend;",file=file,append=TRUE);
    	}
}
