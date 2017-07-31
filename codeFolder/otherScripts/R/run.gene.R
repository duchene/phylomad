run.gene <- function(sdata, format = "phylip", aadata = F, model = "GTR+G", phymlPath, Nsims = 100, para = F, ncore = 1, testStats = c("chisq", "multlik", "delta", "biochemdiv", "consind", "brsup", "CIbrsup", "trlen", "maha"), returnSimPhylo = F, returnSimDat = F, tree = NULL){
	 
	 # Get test statistics
	 
	 if(format == "phylip"){
                  if(aadata) data <- read.aa(sdata) else data <- read.dna(sdata)
         } else if(format == "fasta"){
                  if(aadata) data <- read.aa(sdata, format = "fasta") else data <- read.dna(sdata, format = "fasta")
         } else if(format == "bin"){
                  data <- sdata
         } else if(format == "nexus"){
	   	  if(aadata) data <- as.AAbin(read.nexus.data(sdata)) else data <- as.DNAbin(read.nexus.data(sdata))
	 }
	 
	 empstats <- get.test.statistics(data, format = "bin", aadata = aadata, geneName = "empirical.alignment.phy", phymlPath = phymlPath, model = model, stats = testStats, tree = tree)
	 system("rm empirical.alignment.phy")

	 # Simulate data sets.

	 l <- ncol(data)
	 sim <- list()
	 for(i in 1:Nsims){
	       if(model == "GTR+G"){
               	      rates = phangorn:::discrete.gamma(empstats$alphaParam, k = 4)
               	      rates <- rates + 0.0001
               	      sim_dat_all <- lapply(rates, function(r) simSeq(empstats$outputTree, l = round(l/4, 0), Q = empstats$gtrMatrix, bf = empstats$piParams, rate = r))
               	      sim[[i]] <- as.DNAbin(c(sim_dat_all[[1]], sim_dat_all[[2]], sim_dat_all[[3]], sim_dat_all[[4]]))
	       } else if(model == "GTR"){
	       	      sim[[i]] <- as.DNAbin(simSeq(empstats$outputTree, l = l, Q = empstats$gtrMatrix, bf = empstats$piParams))
	       } else if(model == "HKY+G"){
	       	      rates = phangorn:::discrete.gamma(empstats$alphaParam, k = 4)
                      rates <- rates + 0.0001
                      sim_dat_all <- lapply(rates, function(r) simSeq(empstats$outputTree, l = round(l/4, 0), Q = c(1, 2*empstats$trtvRatio, 1, 1, 2*empstats$trtvRatio, 1), bf = empstats$piParams, rate = r))
                      sim[[i]] <- as.DNAbin(c(sim_dat_all[[1]], sim_dat_all[[2]], sim_dat_all[[3]], sim_dat_all[[4]]))
	       } else if(model == "HKY"){
	       	      sim[[i]] <- as.DNAbin(simSeq(empstats$outputTree, l = l, Q = c(1, 2*empstats$trtvRatio, 1, 1, 2*empstats$trtvRatio, 1), bf = empstats$piParams))
	       } else if(model == "JC+G"){
	       	      rates = phangorn:::discrete.gamma(empstats$alphaParam, k = 4)
                      rates <- rates + 0.0001
		      sim_dat_all <- lapply(rates, function(r) simSeq(empstats$outputTree, l = round(l/4, 0), rate = r))
                      sim[[i]] <- as.DNAbin(c(sim_dat_all[[1]], sim_dat_all[[2]], sim_dat_all[[3]], sim_dat_all[[4]]))
	       } else if(model == "JC"){
	       	      sim[[i]] <- as.DNAbin(simSeq(empstats$outputTree, l = l))
	       } else if(model == "LG+G" | model == "WAG+G" | model == "JTT+G" | model == "Dayhoff+G"){
	       	      simpmod <- gsub("[+]G", "", model)
                      rates = phangorn:::discrete.gamma(empstats$alphaParam, k = 4)
                      rates <- rates + 0.0001
                      sim_dat_all <- lapply(rates, function(r) simSeq(empstats$outputTree, type = "AA", model = simpmod, l = round(l/4, 0), rate = r))
                      sim[[i]] <- as.AAbin(c(sim_dat_all[[1]], sim_dat_all[[2]], sim_dat_all[[3]], sim_dat_all[[4]]))
               } else if(model == "LG" | model == "WAG" | model == "JTT" | model == "Dayhoff"){
                      sim[[i]] <- as.AAbin(simSeq(empstats$outputTree, l = l, type = "AA", model = model))
               }
	 
	 }
	 

	 # Get test statistics for simulations.

	 if(!para){
	   sim.stats <- list()
	 
	   for(i in 1:Nsims){	       
	       sim.stats[[i]] <- get.test.statistics(sim[[i]], format = "bin", aadata = aadata, geneName = paste0("sim.data.", i), phymlPath = phymlPath, model = model, stats = testStats, tree = tree)
	       system(paste0("rm ", paste0("sim.data.", i)))
	   }
	   
	 } else {
	   ### START PARALLEL COMPUTING
	   print('Parallel computing started')
	   require(foreach)
	   require(doParallel)
		
	   runSim <- function(i){
	     tRep <- get.test.statistics(sim[[i]], format = "bin", aadata = aadata, geneName = paste0("sim.data.", i), phymlPath = phymlPath, model = model, stats = testStats, tree = tree)
             system(paste0("rm ", paste0("sim.data.", i)))
	     return(tRep)		
	   }	  
	   cl <- makeCluster(ncore)
	   registerDoParallel(cl)
	   simReps <- foreach(x = 1:Nsims, .packages = c('phangorn', 'ape'), .export = c('get.test.statistics', 'runPhyML', 'get.chisqstat', 'get.biodivstat')) %dopar% runSim(x)
	   sim.stats <- simReps 
	   stopCluster(cl)
	   print("Parallel computing ended successfully")
	  }
	  ### END PARALLEL COMPUTING

	 # Get P-values for test statistics.
	 
	 results <- list()	 

	 if("chisq" %in% testStats){
	 results$emp.chisq <- empstats$chisq
	 results$sim.chisqs <- sapply(sim.stats, function(x) x$chisq)
	 results$chisq.tailp <- length(which(results$sim.chisqs < empstats$chisq)) / Nsims
	 results$chisq.sdpd <- (results$emp.chisq - mean(results$sim.chisqs)) / sd(results$sim.chisqs)
	 }

	 if("multlik" %in% testStats){
	 results$emp.multlik <- empstats$multlik
	 results$sim.multliks <- sapply(sim.stats, function(x) x$multlik)
	 results$multlik.tailp <- length(which(results$sim.multliks < empstats$multlik)) / Nsims
	 results$multlik.sdpd <- (results$emp.multlik - mean(results$sim.multliks)) / sd(results$sim.multliks)
	 }
	 
	 if("delta" %in% testStats){
	 results$emp.delta <- empstats$delta
	 results$sim.deltas <- sapply(sim.stats, function(x) x$delta)
	 results$delta.tailp <- length(which(results$sim.deltas < empstats$delta)) / Nsims
	 results$delta.sdpd <- (results$emp.delta - mean(results$sim.deltas)) /	sd(results$sim.deltas)
	 }
	 
	 if("biochemdiv" %in% testStats){
	 results$emp.biochemdiv <- empstats$biocp
	 results$sim.biochemdiv <- sapply(sim.stats, function(x) x$biocp)
	 results$biochemdiv.tailp <- length(which(results$sim.biochemdiv < empstats$biocp)) / Nsims
	 results$biochemdiv.sdpd <- (results$emp.biochemdiv - mean(results$sim.biochemdiv)) / sd(results$sim.biochemdiv)
	 }
	 
	 if("consind" %in% testStats){
	 results$emp.consind <- empstats$consind
	 results$sim.consind <- sapply(sim.stats, function(x) x$consind)
         results$consind.tailp <- length(which(results$consind < empstats$consind)) / Nsims
	 results$consind.sdpd <- (results$emp.consind - mean(results$sim.consind)) / sd(results$sim.consind)
	 }
	 
	 if("brsup" %in% testStats){
	 results$emp.brsup <- empstats$brsup
	 results$sim.meanbrsup <- sapply(sim.stats, function(x) x$brsup)
	 results$brsup.tailp <- length(which(results$sim.meanbrsup < empstats$brsup)) / Nsims
	 results$brsup.sdpd <- (results$emp.brsup - mean(results$sim.meanbrsup)) / sd(results$sim.meanbrsup)
	 }
	 
	 if("CIbrsup" %in% testStats){
	 results$emp.CIbrsup <- empstats$CIbrsup
	 results$sim.CIbrsup <- sapply(sim.stats, function(x) x$CIbrsup)
	 results$CIbrsup.tailp <- length(which(results$sim.CIbrsup < empstats$CIbrsup)) / Nsims
	 results$CIbrsup.sdpd <- (results$emp.CIbrsup - mean(results$sim.CIbrsup)) / sd(results$sim.CIbrsup)
	 }
	 
	 if("trlen" %in% testStats){
	 results$emp.trlen <- empstats$trlen
	 results$sim.trlens <- sapply(sim.stats, function(x) x$trlen)
	 results$trlen.tailp <- length(which(results$sim.trlens < empstats$trlen)) / Nsims
	 results$trlen.sdpd <- (results$emp.trlen - mean(results$sim.trlens)) /	sd(results$sim.trlens)
	 }

	 if("maha" %in% testStats){
	 all.emp.stats <- unlist(results[grep("emp[.]", names(results))])
	 all.sim.stats <- do.call(cbind, results[grep("sim[.]", names(results))])
	 all.stats.mat <- rbind(all.sim.stats, all.emp.stats)
	 failstats <- vector()
	 if(length(all.emp.stats) == 0){
	 print("No test statistics were calculated.")
	 } else {
	 for(i in 1:ncol(all.stats.mat)){
	        if(length(unique(all.sim.stats[,i])) == 1 || any(is.na(all.sim.stats[,i]))){
			print(paste(colnames(all.stats.mat)[i], "was not included in the mahalanobis calculation and is likely to be unreliable. This is possibly because the values for all simulations are the same."))
	 		failstats <- c(failstats, i)
	 	}
	 }
	 if(length(failstats) > 0){
	 	all.sim.stats <- all.sim.stats[,-failstats]
           	all.stats.mat <- all.stats.mat[,-failstats]
	 }
	 #print(all.stats.mat)
	 if(length(all.stats.mat) >= 2*Nsims){
	 	mahavector <- mahalanobis(all.stats.mat, colMeans(all.stats.mat, na.rm = T), cov(all.stats.mat, use = "complete.obs"))
	 	results$emp.maha <- tail(mahavector, 1)
	 	results$sim.maha <- mahavector[1:Nsims]
	 	results$maha.tailp <- length(which(mahavector[1:Nsims] > results$emp.maha)) / Nsims
	 	results$maha.sdpd <- (results$emp.maha - mean(results$sim.maha)) / sd(results$sim.maha)
	 } else {
	   	print("Mahalanobis cannot be calculated. This is possibly because all the statistics produced results that cannot be used.")
	 }
	 }
	 }
	 
	 results$empirical.tree <- empstats$outputTree

	 if(length(grep("HKY|GTR", model)) == 1) results$piParams <- empstats$piParams
	 
	 if(length(grep("[+]G", model)) == 1) results$alphaParam <- empstats$alphaParam
	 
	 if(length(grep("GTR", model)) == 1) results$gtrMatrix <- empstats$gtrMatrix
	 
         if(length(grep("HKY", model)) == 1) results$trtvRatio <- empstats$trtvRatio
	 
	 if(returnSimPhylo){
		results$simPhylos <- lapply(sim.stats, function(x) x$outputTree)
		class(results$simPhylos) <- "multiPhylo"
	 }
	 if(returnSimDat) results$simDat <- sim

	 return(results)

}