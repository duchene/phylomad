`genetree.vector` <-
function(filenames,outputfile)
{
    nfile<-length(filenames)

    treestr<-read.tree.string(filenames[1])$tree
    for(i in 2:nfile){
	treestr1<-read.tree.string(filenames[i])$tree
	treestr<-cbind(treestr,treestr1)
    }

    treestr1<-as.vector(treestr)
    write.tree.string(treestr1,outputfile,format="phylip")
    return(treestr)

}
