
test.saturation <- function(sdata = loci, format = "phylip", para = parallelise, ncore = 1, clean = "cleandata", plotdat = F){

	 runLoc <- function(loc, format = format, aadata = aadata, clean = clean){
             if(clean == "cleandata") cleanlogical <- T else cleanlogical <- F
	     genebin <- clean.gene(sdata = loc, format = format, aadata = F, clean = cleanlogical, plotdat = plotdat)
	     
	     ## FIGURE OUT CODON POSITION TEST
	     if(clean == "codonpos")
	     
	     ## FILL IN CODE TO TEST SATURATION!!

	     satres <- "banana"
	     
	     ##
	     
	     dist.raw <- dist.dna(loc, model = "raw", pairwise.deletion = T)
	     dist.tn93 <- dist.dna(loc, model = "TN93", pairwise.deletion = T)
	     distdat <- list(dist.raw, dist.tn93)
	     
             tRep <- list(saturation.test.results = satres, genetic.distance.data = distdat)
             return(tRep)
         }

	 reslist <- list()

         if(!para){

           for(i in 1:length(loci)){
	       reslist[[i]] <- runLoc(loci[i], format = format, aadata = F, clean = clean, plotdat = plotdat)
           }

         } else {
           ### START PARALLEL COMPUTING
           print('Parallel computing started')
           require(foreach)
           require(doParallel)
           cl <- makeCluster(ncore)
           registerDoParallel(cl)
           reslist <- foreach(x = 1:length(loci), .packages = c('phangorn', 'ape'), .export = c('format', 'aadata', 'clean', 'clean.gene')) %dopar% runSim(x, format = format, aadata = F, clean = clean, plotdat = plotdat)
           stopCluster(cl)
           print("Parallel computing ended successfully")
          }
          ### END PARALLEL COMPUTING
	  
	  names(reslist) <- loci

	  return(reslist)

}