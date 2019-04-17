`species.name` <-
function(str)
{
   	str<-gsub("[ ]", "", str)
	str<-gsub("[;.]","",str)
	str<-gsub(":[0-9e-]*","",str)
	str<-gsub("#[0-9e-]*","",str)
	str<-gsub(")[0-9e-]*","",str)
	str<-gsub("[()]","",str)
	str<-gsub("[ ]", "", str)
	name<-sort(unlist(strsplit(str,split=",")))

    	return(name)
}

