`sim.coaltree.sp` <-
function(rootnode,nodematrix,nspecies,seq,name)
{
    theta<-nodematrix[rootnode,5]

    if(rootnode<=nspecies){
	{if(seq[rootnode] == 1){
		z<-list(gt="", height=as.matrix, node=as.matrix)
		z$gt<-name[rootnode]
		z$height<-0
		z$node<-nodematrix
		return(z)}
	 else{
		treestr<-paste(name[rootnode],"s",1:seq[rootnode],sep="")
		i<-seq[rootnode]
		height<-rexp(1,rate=i*(i-1)/theta)
    		brlens<-rep(0,i)
		father<-nodematrix[rootnode,1]
		fatherheight<-node.height(father,nodematrix,nspecies)
		
		while(height<fatherheight){
			nodematrix[rootnode,6]<-nodematrix[rootnode,6]+1
			##randomly choose two nodes
 			b<-sample(1:i,2)
		
			##update groups
    			newname<-paste("(",treestr[b[1]],sep="")
    			newname<-paste(newname,":",sep="")
    			newname<-paste(newname,round(height-brlens[b[1]],6),sep="")
    			newname<-paste(newname,",",sep="")   
    			newname<-paste(newname,treestr[b[2]],sep="")
    			newname<-paste(newname,":",sep="")
    			newname<-paste(newname,round(height-brlens[b[2]],6),sep="")
    			newname<-paste(newname,")",sep="")
    			treestr[b[1]]<-newname
			brlens[b[1]]<-height

			##update dist,treestr, and branch length
			index<-1:i
    			index[b[2]]<-0
    			index<-index[index>0]
        		treestr<-treestr[index]
			brlens<-brlens[index]
			if(i==2)
				break
			i<-i-1
			height<-height+rexp(1,rate=i*(i-1)/theta)
		}
		z<-list(gt="", height=as.matrix,node=as.matrix)
		z$gt<-treestr
		z$height<-brlens
		z$node<-nodematrix
		return(z)
	}
       }
    	
    }
    if(rootnode>nspecies){
	son1<-nodematrix[rootnode,2]
	son2<-nodematrix[rootnode,3]
	leftstr<-sim.coaltree.sp(rootnode=son1,nodematrix=nodematrix,nspecies,seq,name)
	nodematrix<-leftstr$node
	rightstr<-sim.coaltree.sp(rootnode=son2,nodematrix=nodematrix,nspecies,seq,name)
	nodematrix<-rightstr$node
	i<-length(leftstr$gt)+length(rightstr$gt)

	treestr<-1:i
	treestr[1:length(leftstr$gt)]<-leftstr$gt
	treestr[(length(leftstr$gt)+1):i]<-rightstr$gt

 	brlens<-1:i
	brlens[1:length(leftstr$height)]<-leftstr$height
	brlens[(length(leftstr$height)+1):i]<-rightstr$height



	height<-rexp(1,rate=i*(i-1)/theta)+node.height(rootnode,nodematrix,nspecies)
	father<-nodematrix[rootnode,1]
	if(father == -9 | father == -8)
		fatherheight<-100000 else
		fatherheight<-node.height(father,nodematrix,nspecies)
    	brlens

	
	while(height<fatherheight){
		nodematrix[rootnode,6]<-nodematrix[rootnode,6]+1
		##randomly choose two nodes
 		b<-sample(1:i,2)
		
		##update groups
    		newname<-paste("(",treestr[b[1]],sep="")
    		newname<-paste(newname,":",sep="")
    		newname<-paste(newname,round(height-brlens[b[1]],6),sep="")
    		newname<-paste(newname,",",sep="")   
    		newname<-paste(newname,treestr[b[2]],sep="")
    		newname<-paste(newname,":",sep="")
    		newname<-paste(newname,round(height-brlens[b[2]],6),sep="")
    		newname<-paste(newname,")",sep="")
    		treestr[b[1]]<-newname
		brlens[b[1]]<-height

		##update dist,treestr, and branch length
		index<-1:i
    		index[b[2]]<-0
    		index<-index[index>0]
        	treestr<-treestr[index]
		brlens<-brlens[index]
		if(i==2)
			break
		i<-i-1
		height<-height+rexp(1,rate=i*(i-1)/theta)
	}
	if(nodematrix[rootnode,1]==-9 | nodematrix[rootnode,1]==-8)
		treestr<-paste(treestr,";",sep="")
	z<-list(gt="", height=as.matrix,node=as.matrix)
	z$gt<-treestr
	z$node<-nodematrix
	z$height<-brlens
    	return(z)
    }

}