
test.saturation <- function(loci, format = "phylip", para = parallelise, ncore = 1, clean = "cleandata", stats = stats, plotdat = F){
	 aadata <- F
	 runLoc <- function(loc, format = format, aadata = aadata, clean = clean, stats = stats, plotdat = plotdat){
             if(clean == "cleandata") cleanlogical <- T else cleanlogical <- F
	     genebin <- list(clean.gene(sdata = loc, format = format, aadata = aadata, clean = cleanlogical))
	     
	     if(clean == "codonpos") genebin <- list(pos1and2 = genebin[[1]][,c(T,T,F)], pos3 = genebin[[1]][,c(F,F,T)])
	     
	     if(plotdat) distdat <- list()
	     
	     genesres <- list()
	     
	     for(i in 1:length(genebin)){
	     
		#### TESTS OF SATURATION ####
		
		satres <- list() 
        	colnames(saturation.result[[i]]) <- c("t", "t_moderate_risk_threshold", "t_high_risk_threshold", "Risk")
		
		if("cith" %in% stats | "comth" %in% stats) phymlres <- runPhyML(genebin[[i]], format = format, aadata = aadata, temp_name = "empirical", phymlPath = phymlPath, model = "GTR+G")

        	if("enth" %in% stats){
                	  satres$enres <- vector()
			  satres$enres[1] <- get.entropy.test(al)$t
			  satres$enres[2] <- getSatThreshold(ncol(genebin[[i]]), nrow(genebin[[i]]), th1entsqrt)
			  satres$enres[3] <- getSatThreshold(ncol(genebin[[i]]), nrow(genebin[[i]]), th2entsqrt)
		}
        	if("cith" %in% stats){
			  satres$cires <- vector()
			  satres$cires[1] <- get.ci.test(al, phymlres$tree)$t
			  satres$cires[2] <- getSatThreshold(ncol(genebin[[i]]), nrow(genebin[[i]]), th1entsqrt)
                          satres$cires[3] <- getSatThreshold(ncol(genebin[[i]]), nrow(genebin[[i]]), th2entsqrt)
        	}
        	if("comth" %in% stats){
			   satres$comres <- vector()
                	   satres$comres[1] <- get.comp.test(al, phymlres$tree)$t
			   satres$comres[2] <- getSatThreshold(ncol(genebin[[i]]), nrow(genebin[[i]]), th1entsqrt)
                           satres$comres[3] <- getSatThreshold(ncol(genebin[[i]]), nrow(genebin[[i]]), th2entsqrt)
        	}
		
		for(j in 1:length(satres)) if(satres[[j]][1] > satres[[j]][2]) satres[[j]][4] <- "LOW" else if(satres[[j]][1] < satres[[j]][3]) satres[[j]][4] <- "HIGH" else satres[[j]][4] <- "MEDIUM"

		genesres[[i]] <- do.call(c, satres)
		names(genesres[[i]]) <- as.character(sapply(stats, function(x) paste(c("t", "t_moderate_risk_threshold", "t_high_risk_threshold", "Risk"), x, sep = "_")))
	     	
		###########################################
	     	
		if(plotdat){
			dist.raw <- dist.dna(genebin[[i]], model = "raw", pairwise.deletion = T)
	     		dist.tn93 <- dist.dna(genebin[[i]], model = "TN93", pairwise.deletion = T)
	     		distdat[[i]] <- list(dist.raw, dist.tn93)
	     	}
	     }
	     
	     if(length(genesres) > 1) genesrestab <- do.call(rbind, genesres) else genesrestab <- matrix(genesres[[1]], 1, 4)
	     
	     if(length(genebin) > 1) rownames(genesrestab) <- paste(loc, names(genebin), sep = "_") else rownames(genesrestab) <- loc
	     
             if(plotdat) tRep <- list(saturation.test.results = genesres, genetic.distance.data = distdat) else tRep <- list(saturation.test.results = satres)
             return(tRep)
         }

	 reslist <- list()

         if(!para){

           for(i in 1:length(loci)){
	       reslist[[i]] <- runLoc(loci[i], format = format, aadata = aadata, clean = clean, stats = stats, plotdat = plotdat)
           }

         } else {

           ### START PARALLEL COMPUTING
           print('Parallel computing started')
           require(foreach)
           require(doParallel)
           cl <- makeCluster(ncore)
           registerDoParallel(cl)
           reslist <- foreach(x = loci, .packages = c('phangorn', 'ape'), .export = c('clean.gene')) %dopar% runLoc(x, format = format, aadata = aadata, clean = clean, stats = stats, plotdat = plotdat)
           stopCluster(cl)
           print("Parallel computing ended successfully")
           ### END PARALLEL COMPUTING

	  }
	  
	  if(length(reslist) > 1) restab <- do.call(rbind, unname(reslist)) else restab <- reslist[[1]]
	  
	  return(restab)
}