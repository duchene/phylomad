
test.saturation <- function(loci, format = "phylip", iqtreePath = iqtreePath, para = parallelise, ncore = 1, clean = "cleandata", stats = stats, plotdat = F, linmods = funclist){
	 
	 aadata <- F
	 
	 runLoc <- function(loc, format = format, aadata = aadata, clean = clean, stats = stats, iqtreePath = iqtreePath, plotdat = plotdat, satfunclist = linmods){
             if(clean == "cleandata") cleanlogical <- T else cleanlogical <- F
	     genebin <- list(clean.gene(sdata = loc, format = format, aadata = aadata, clean = cleanlogical))
	     
	     if(clean == "codonpos") genebin <- list(pos1and2 = genebin[[1]][,c(T,T,F)], pos3 = genebin[[1]][,c(F,F,T)])
	     
	     if(plotdat) distdat <- list()
	     
	     genesres <- list()
	     
	     getSatThreshold <- function(seqlen, ntax, satthresfun){
		 resthres <- predict(satthresfun, newdata = data.frame(seqlen = seqlen, ntax = ntax))
		 return(resthres)
	     }
	     
	     for(i in 1:length(genebin)){
	     
		if("cith" %in% stats | "comth" %in% stats){
			dist.tn93 <- dist.dna(genebin[[i]], model = "TN93", pairwise.deletion = T)
			tree <- njs(dist.tn93)
		}
		
                if(plotdat){
                        dist.raw <- dist.dna(genebin[[i]], model = "raw", pairwise.deletion = T)
                        if(!"cith" %in% stats | !"comth" %in% stats) dist.tn93 <- dist.dna(genebin[[i]], model = "TN93", pairwise.deletion = T)
                        distdat[[i]] <- list(dist.raw, dist.tn93)
		}

		#### TESTS OF SATURATION ####
		
		satres <- list()
		
#		if("cith" %in% stats | "comth" %in% stats){
#			  tempalname <- paste(sample(letters, 5), collapse = "")
#			  iqtreeres <- runIQtree(genebin[[i]], format = "bin", aadata = aadata, temp_name = tempalname, iqtreePath = iqtreePath, model = "GTR+G")
#			  system(paste("rm", tempalname))
#			  tree <- iqtreeres$tree
#		}
		
		Nsites <- ncol(genebin[[i]])
		Ntax <- nrow(genebin[[i]])
		
        	if("enth" %in% stats){
                	  satres$enres <- vector()
			  satres$enres[1] <- get.entropy.test(genebin[[i]])$t
			  satres$enres[2] <- getSatThreshold(Nsites, Ntax, satfunclist[["th1entsqrt"]])
			  satres$enres[3] <- getSatThreshold(Nsites, Ntax, satfunclist[["th2entsqrt"]])
			  satres$enres <- round(satres$enres, 2)
		}
        	if("cith" %in% stats){
			  satres$cires <- vector()
			  satres$cires[1] <- get.ci.test(tree, genebin[[i]])$t * (-1)
			  satres$cires[2] <- getSatThreshold(Nsites, Ntax, satfunclist[["th1cisqrt"]])
                          satres$cires[3] <- getSatThreshold(Nsites, Ntax, satfunclist[["th2cisqrt"]])
			  satres$cires <- round(satres$cires, 2)
        	}
        	if("comth" %in% stats){
			   satres$comres <- vector()
                	   satres$comres[1] <- get.comp.test(tree, genebin[[i]])$t * (-1)
			   satres$comres[2] <- getSatThreshold(Nsites, Ntax, satfunclist[["th1comsqrt"]])
                           satres$comres[3] <- getSatThreshold(Nsites, Ntax, satfunclist[["th2comsqrt"]])
			   satres$comres <- round(satres$comres, 2)
        	}
		
		for(j in 1:length(satres)) if(satres[[j]][1] > satres[[j]][2]) satres[[j]][4] <- "LOW" else if(satres[[j]][1] < satres[[j]][3]) satres[[j]][4] <- "HIGH" else satres[[j]][4] <- "MEDIUM"

		genesres[[i]] <- c(Nsites, Ntax, do.call(c, satres))
		names(genesres[[i]]) <- c("N_sites", "N_taxa", as.character(sapply(stats, function(x) paste(c("t", "t_moderate_risk_threshold", "t_high_risk_threshold", "Risk"), x, sep = "_"))))
		
		###########################################

	     }
	     
	     if(length(genesres) > 1){
	     		genesrestab <- do.call(rbind, genesres)
	     } else {
			genesrestab <- matrix(genesres[[1]], 1, 2+(length(satres)*4))
			colnames(genesrestab) <- names(genesres[[1]])
	     }
	     
             if(plotdat) tRep <- list(saturation.test.results = genesrestab, genetic.distance.data = distdat) else tRep <- list(saturation.test.results = genesrestab)
	     
             return(tRep)
         }

	 reslist <- list()

         if(!para){

           for(i in 1:length(loci)){
	       reslist[[i]] <- runLoc(loci[i], format = format, aadata = aadata, clean = clean, stats = stats, iqtreePath = iqtreePath, plotdat = plotdat, satfunclist = linmods)
           }

         } else {

           ### START PARALLEL COMPUTING
           print('Parallel computing started')
           require(foreach)
           require(doParallel)
           cl <- makeCluster(ncore)
           registerDoParallel(cl)
           reslist <- foreach(x = loci, .packages = c('phangorn', 'ape'), .export = c('clean.gene', 'linmods', 'runIQtree', 'get.entropy.test', 'get.ci.test', 'get.comp.test')) %dopar% runLoc(x, format = format, aadata = aadata, clean = clean, stats = stats, iqtreePath = iqtreePath, plotdat = plotdat, satfunclist = linmods)
           stopCluster(cl)
           print("Parallel computing ended successfully")
           ### END PARALLEL COMPUTING

	  }
	  
	  return(reslist)
}