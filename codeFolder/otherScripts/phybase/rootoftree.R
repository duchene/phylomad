`rootoftree` <-
function (nodematrix)
{
    if(sum(nodematrix[,1]<0)>1)
	stop("more than two roots!")
    nodes<-sum((nodematrix[,1]<0)*(1:length(nodematrix[,1])))
    return(nodes)
}

