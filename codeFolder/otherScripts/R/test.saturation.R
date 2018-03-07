
test.saturation <- function(loci, format = "phylip", para = parallelise, ncore = 1, clean = "cleandata", plotdat = F){
	 aadata <- F
	 runLoc <- function(loc, format = format, aadata = aadata, clean = clean, plotdat = plotdat){
             if(clean == "cleandata") cleanlogical <- T else cleanlogical <- F
	     genebin <- list(clean.gene(sdata = loc, format = format, aadata = F, clean = cleanlogical))
	     
	     if(clean == "codonpos") genebin <- list(pos1and2 = genebin[[1]][,c(T,T,F)], pos3 = genebin[[1]][,c(F,F,T)])
	     
	     satres <- list()
	     if(plotdat) distdat <- list()
	     
	     for(i in 1:length(genebin)){
	     
		## FILL IN CODE TO TEST SATURATION!!

	     	satres[[i]] <- "banana"
	     
		##
	     	
		if(plotdat){
			dist.raw <- dist.dna(genebin[[i]], model = "raw", pairwise.deletion = T)
	     		dist.tn93 <- dist.dna(genebin[[i]], model = "TN93", pairwise.deletion = T)
	     		distdat[[i]] <- list(dist.raw, dist.tn93)
	     	}
	     }
	     
             if(plotdat) tRep <- list(saturation.test.results = satres, genetic.distance.data = distdat) else tRep <- list(saturation.test.results = satres)
             return(tRep)
         }

	 reslist <- list()

         if(!para){

           for(i in 1:length(loci)){
	       reslist[[i]] <- runLoc(loci[i], format = format, aadata = aadata, clean = clean, plotdat = plotdat)
           }

         } else {

           ### START PARALLEL COMPUTING
           print('Parallel computing started')
           require(foreach)
           require(doParallel)
           cl <- makeCluster(ncore)
           registerDoParallel(cl)
           reslist <- foreach(x = loci, .packages = c('phangorn', 'ape'), .export = c('clean.gene')) %dopar% runLoc(x, format = format, aadata = aadata, clean = clean, plotdat = plotdat)
           stopCluster(cl)
           print("Parallel computing ended successfully")
           ### END PARALLEL COMPUTING

	  }
	  
	  names(reslist) <- loci
	  return(reslist)
}