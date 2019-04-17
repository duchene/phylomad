`node2name` <-
function(treestr, name="")
{
    str<-treestr    
    speciesname<-species.name(str)
    

    old<-paste(1:length(name),":",sep="")
    new<-paste(name,"$",sep="")
    for(i in length(name):1)	
	str<-gsub(old[i],new[i],str)
    str<-gsub("\\$",":",str)

    str

}