`write.tree.string` <-
function (X,format="Nexus",file="",name="")
{

   if(length(name)==1)
	speciesname<-species.name(X[1])
    else
	speciesname<-name

    if(length(grep("nexus", format, ignore.case = TRUE))>0){
	write("#NEXUS",file)
	write("begin trees;",file,append=TRUE)
	write("  translate",file,append=TRUE)
	for(i in 1:(length(speciesname)-1)){
		a<-paste(i,speciesname[i])
		a<-paste(a,",",sep="")
		a<-paste("\t",a)
		write(a,file,append=TRUE)
	}
	a<-paste(length(speciesname),speciesname[length(speciesname)])
	a<-paste(a,";",sep="")
	a<-paste("\t",a)
	write(a,file,append=TRUE)
    	for(i in 1:length(X)){
		a<-paste("  tree",i)
		a<-paste(a,"=")
		X[i]<-name2node(X[i],speciesname)
		a<-paste(a,X[i])
		write(a,file,append=TRUE)
    	}
    	write("end;",file,append=TRUE)
    }


    if(length(grep("phylip", format, ignore.case = TRUE))>0)
 	write(X,file)

}

