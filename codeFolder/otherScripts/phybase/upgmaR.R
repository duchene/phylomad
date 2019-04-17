`upgmaR` <-
 
 # This function takes a distance matrix and produces an ultrametric tree using UPGMA.
 # NOTE: method "min" from phybase has been REMOVED
 
  function(dist,name,method="average")
  {
    nspecies<-length(name)       # Number of species
    treestr<-name                       # Names to be used and modified to build Newick format tree.
    brlens<-rep(0,nspecies)        # Vector of distance from leaves to MRCA of clade (initially zero for 1-clades)
    cladesize<- rep(1, nspecies) # Vector of counts of number of taxa in each clade already grouped (initially all 1s). 
    b<-rep(0,2)                            # This will be location of the smallest positive entry in dist 
      
   
    for(i in 1:(nspecies-1)){
        ## Find the minimum distance, determining clades to be amalgamated
        diag(dist)<-NA                                                               # Change the diagonal of the distance matrix to NA so that the minimum number is not 0
        position <- which(dist == min(dist, na.rm = TRUE))[1]  # Find the the first occurrence of the smallest number, searching column-wise.
        b[1] <- floor(position/length(dist[1, ]))                            # This will be the column index.
        b[2] <- position - b[1] * length(dist[1, ])                          # This will be the row  index.
        if (b[2] == 0)                                                                  # If the position is at the bottom of the column, then we must adjust
            b[2] <- length(dist[1, ]) 
        else b[1] <- b[1] + 1
     
       ## Amalgamate clades and collapse to smaller distance matrix
       halfdist12= dist[ b[2], b[1] ]/2    # Save half the smallest value in the matrix
      
       if(method == "average"){    
           dist[b[1],]= ( dist[b[1],]*cladesize[b[1]] + dist[b[2],]*cladesize[b[2]])/(cladesize[b[1]] +cladesize[b[2]]) # Compute weighted average using clade size for new row/column of distance matrix
           dist[,b[1]]= t(dist[b[1],])                                                                                                                         
       }
       if (method == "min") {
            bdist <- dist[b[2], b[1]]
            dist[, b[1]] <- pmin(dist[, b[1]], dist[, b[2]])
            dist[b[1], ] <- pmin(dist[b[1], ], dist[b[2], ])          
       }
       # Remove a column and row
       index <- 1:(nspecies + 1 - i)
       index[b[2]]<-0
       index<-index[index>0]
       dist<-dist[index,index]
                
      ##Build a Newick tree
      newname<-paste("(",treestr[b[1]],sep="")
      newname<-paste(newname,":",sep="")
      newname<-paste(newname,round((halfdist12-brlens[b[1]]),5),sep="")
      newname<-paste(newname,",",sep="")
      newname<-paste(newname,treestr[b[2]],sep="")
      newname<-paste(newname,":",sep="")
      newname<-paste(newname,round((halfdist12-brlens[b[2]]),5),sep="")
      newname<-paste(newname,")",sep="")
      treestr[b[1]]<-newname
      brlens[b]<-halfdist12
           
      ##Update clade sizes.
      cladesize[b[1]]<- cladesize[b[1]] + cladesize[b[2]]
      cladesize[b[2]]<- 0
      cladesize<- cladesize[cladesize>0]
     
      ##
      treestr<-treestr[index]
      brlens<-brlens[index]
     }
   
    ## Display results
    treestr<-paste(treestr,";",sep="")
    z<-list(nodes=matrix,treestr="",name="")
    node<-read.tree.nodes(treestr,name)
    z$nodes<-node$nodes
    z$name<-node$name
    z$treestr<-treestr
    z
  }

