`noclock2clock` <-
function(inode, treematrix, nspecies)
{
	if(inode > nspecies)
	{
		son1<-treematrix[inode,2]
   		son2<-treematrix[inode,3]
		treematrix<-noclock2clock(son1,treematrix,nspecies)
		treematrix<-noclock2clock(son2,treematrix,nspecies)

    		leftheight<-treematrix[son1,4]+node.height(son1,treematrix,nspecies)
		rightheight<-treematrix[son2,4]+node.height(son2,treematrix,nspecies)
 		leftratio<-(leftheight+rightheight)/2/leftheight
		rightratio<-(leftheight+rightheight)/2/rightheight

   		leftsonnodes<-offspring.nodes(son1,treematrix,nspecies)
   		rightsonnodes<-offspring.nodes(son2,treematrix,nspecies)

        	treematrix[leftsonnodes,4]<-treematrix[leftsonnodes,4]*leftratio
		treematrix[rightsonnodes,4]<-treematrix[rightsonnodes,4]*rightratio
    	}
    	return(treematrix)	
}