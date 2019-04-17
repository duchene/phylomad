'bootstrap.mulgene' <-
function(sequence,gene,name,boot,outfile=""){
    write.table("",file=outfile,col.names=FALSE,row.names=FALSE,quote=FALSE)
    #row:taxa, col:sites
    sequence<-t(sequence)
    ntaxa<-dim(sequence)[2]
    for(j in 1:boot){
	ngene<-dim(gene)[1]
	index<-sample(1:ngene,ngene,replace=TRUE)
	bootgene<-gene[index,]
	for(i in 1:ngene){
		gene1<-sequence[bootgene[i,1]:bootgene[i,2],]
		index<-sample(1:dim(gene1)[1],dim(gene1)[1],replace=TRUE)
		gene2<-gene1[index,]
		
		#delete the missing taxa
		missing<-rep(1,ntaxa)
		for(k in 1:ntaxa)
		{
		if((sum(gene2[,k]=="?")+sum(gene2[,k]=="-")+sum(gene2[,k]=="N")+sum(gene2[,k]=="n")) == length(gene2[,k])) missing[k]<-0
		}
		nomissing<-(1:ntaxa)*missing
		nomissing<-nomissing[nomissing>0]
		write.dna(t(gene2[,nomissing]),name[nomissing],format="phylip",file=outfile,append=TRUE)
	}
    }
}