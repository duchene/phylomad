test.saturation <- function(loci, iqtreePath = iqtreePath, para = parallelise, ncore = 1, clean = "cleandata", stats = stats, plotdat = F, linmods = funclist){
	 
	 aadata <- F
	 
	 runLoc <- function(loc, aadata = aadata, clean = clean, stats = stats, iqtreePath = iqtreePath, plotdat = plotdat, satfunclist = linmods){
             if(clean == "cleandata") cleanlogical <- T else cleanlogical <- F
	     
	     firstLine <- readLines(loc, n = 1)
if(grepl("[>]", firstLine)) dataFormat <- "fasta" else if(grepl("[#]NEXUS|[#]nexus", firstLine)) dataFormat <- "nexus" else dataFormat <- "phylip"
	     
	     genebin <- list(clean.gene(sdata = loc, format = dataFormat, aadata = aadata, clean = cleanlogical))
	     
	     if(clean == "codonpos") genebin <- list(pos1and2 = genebin[[1]][,c(T,T,F)], pos3 = genebin[[1]][,c(F,F,T)])
	     
	     if(plotdat) distdat <- list()
	     
	     genesres <- list()
	     
	     getSatThreshold <- function(seqlen, ntax, satthresfun){
		 resthres <- predict(satthresfun, newdata = data.frame(seqlen = seqlen, ntax = ntax))
		 return(resthres)
	     }
	     
	     for(i in 1:length(genebin)){
		
                if(plotdat){
                        dist.raw <- dist.dna(genebin[[i]], model = "raw", pairwise.deletion = T)
              		dist.tn93 <- dist.dna(genebin[[i]], model = "TN93", pairwise.deletion = T)
                        distdat[[i]] <- list(dist.raw, dist.tn93)
		}

		#### TESTS OF SATURATION ####
		
		satres <- list()
		
		Nsites <- ncol(genebin[[i]])
		Ntax <- nrow(genebin[[i]])
		
        	if("enth" %in% stats){
                	  satres$enres <- vector()
			  satres$enres[1] <- get.entropy.test(genebin[[i]])$t
			  satres$enres[2] <- getSatThreshold(Nsites, Ntax, satfunclist[[1]])
			  satres$enres[3] <- getSatThreshold(Nsites, Ntax, satfunclist[[3]])
			  satres$enres[4] <- getSatThreshold(Nsites, Ntax, satfunclist[[5]])
			  satres$enres <- round(satres$enres, 2)
		}
        	if("enthvar" %in% stats){
			  satres$envarres <- vector()
			  satres$envarres[1] <- get.entropy.test(genebin[[i]], only.varsites = T)$t
			  satres$envarres[2] <- getSatThreshold(Nsites, Ntax, satfunclist[[2]])
                          satres$envarres[3] <- getSatThreshold(Nsites, Ntax, satfunclist[[4]])
			  satres$envarres[4] <- getSatThreshold(Nsites, Ntax, satfunclist[[6]])
			  satres$envarres <- round(satres$envarres, 2)
        	}
		
		for(j in 1:length(satres)) if(satres[[j]][1] > satres[[j]][2]) satres[[j]][5] <- "low.risk" else satres[[j]][5] <- "high.risk"

		genesres[[i]] <- c(Nsites, Ntax, do.call(c, satres))
		names(genesres[[i]]) <- c("N_sites", "N_taxa", as.character(sapply(stats, function(x) paste(c("t", "t_predicted_threshold", "t_predicted_TPR", "t_predicted_FPR", "Risk"), x, sep = "_"))))
		
		###########################################

	     }
	     
	     if(length(genesres) > 1){
	     		genesrestab <- do.call(rbind, genesres)
	     } else {
			genesrestab <- matrix(genesres[[1]], 1, 2+(length(satres)*5))
			colnames(genesrestab) <- names(genesres[[1]])
	     }
	     
             if(plotdat) tRep <- list(saturation.test.results = genesrestab, genetic.distance.data = distdat) else tRep <- list(saturation.test.results = genesrestab)
	     
             return(tRep)
         }

	 reslist <- list()

         if(!para){

           for(i in 1:length(loci)){
	       reslist[[i]] <- try(runLoc(loci[i], aadata = aadata, clean = clean, stats = stats, iqtreePath = iqtreePath, plotdat = plotdat, satfunclist = linmods))
	       if(class(reslist[[i]]) == "try-error") next
           }

         } else {

           ### START PARALLEL COMPUTING
           print('Parallel computing started')
           require(foreach)
           require(doParallel)
           cl <- makeCluster(ncore)
           registerDoParallel(cl)
           reslist <- foreach(x = loci, .packages = c('phangorn', 'ape'), .export = c('clean.gene', 'linmods', 'runIQtree', 'get.entropy.test')) %dopar% tryCatch(runLoc(x, aadata = aadata, clean = clean, stats = stats, iqtreePath = iqtreePath, plotdat = plotdat, satfunclist = linmods), error = function(e) NULL)
           stopCluster(cl)
           print("Parallel computing ended successfully")
           ### END PARALLEL COMPUTING

	  }
	  
	  return(reslist)
}