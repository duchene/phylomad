run.gene <- function(sdata, format = "phylip", model = "GTR+G", phymlPath, Nsims = 100, para = F, ncore = 1, testStats = c("chisq", "multlik", "delta", "biochemdiv", "consind", "brsup", "trlen", "maha")){
	 
	 # Get test statistics
	 
	 if(format == "phylip"){
                  data <- read.dna(sdata)
         } else if(format == "fasta"){
                  data <- read.dna(sdata, format = "fasta")
         } else if(format == "DNAbin"){
                  data <- sdata
         } else if(format == "nexus"){
	   	  data <- as.DNAbin(read.nexus.data(sdata))
	 }
	 
	 empstats <- get.test.statistics(data, format = format, geneName = "empirical", phymlPath = phymlPath, model = model, stats = testStats)

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
	       	      
	       	      ### HKY MODELS NEED TO BE FINISHED
	       	      rates = phangorn:::discrete.gamma(empstats$alphaParam, k = 4)
                      rates <- rates + 0.0001
                      sim_dat_all <- lapply(rates, function(r) simSeq(empstats$outputTree, l = round(l/4, 0), Q = empstats$gtrMatrix, bf = empstats$piParams, rate = r))
                      sim[[i]] <- as.DNAbin(c(sim_dat_all[[1]], sim_dat_all[[2]], sim_dat_all[[3]], sim_dat_all[[4]]))
	       } else if(model == "HKY"){
	       	      sim[[i]] <- as.DNAbin(simSeq(empstats$outputTree, l = l, Q = empstats$gtrMatrix, bf = empstats$piParams))
		      
	       } else if(model == "JC+G"){
	       	      rates = phangorn:::discrete.gamma(empstats$alphaParam, k = 4)
                      rates <- rates + 0.0001
		      sim_dat_all <- lapply(rates, function(r) simSeq(empstats$outputTree, l = round(l/4, 0), rate = r))
                      sim[[i]] <- as.DNAbin(c(sim_dat_all[[1]], sim_dat_all[[2]], sim_dat_all[[3]], sim_dat_all[[4]]))
	       } else if(model == "JC"){
	       	      sim[[i]] <- as.DNAbin(simSeq(empstats$outputTree, l = l))
	       }
	 
	 }
	 

	 # Get test statistics for simulations.

	 if(!para){
	   sim.stats <- list()
	 
	   for(i in 1:Nsims){	       
	       sim.stats[[i]] <- get.test.statistics(sim[[i]], format = "DNAbin", geneName = paste0("sim.data.", i), phymlPath = phymlPath, model = model, stats = testStats)
	       system(paste0("rm ", paste0("sim.data.", i)))
	   }
	   
	 } else {
	   ### START PARALLEL COMPUTING
	   print('Enter parallel computing')
	   require(foreach)
	   require(doParallel)
		
	   runSim <- function(i){
	     tRep <- get.test.statistics(sim[[i]], format = "DNAbin", geneName = paste0("sim.data.", i), phymlPath = phymlPath, model = model, stats = testStats)
             system(paste0("rm ", paste0("sim.data.", i)))
	     return(tRep)		
	   }	  
	   cl <- makeCluster(ncore)
	   registerDoParallel(cl)
	   simReps <- foreach(x = 1:Nsims, .packages = c('phangorn', 'ape'), .export = c('get.test.statistics', 'runPhyML', 'getchisqs')) %dopar% runSim(x)
	   sim.stats <- simReps 
	   stopCluster(cl)
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
	 results$sim.multinoms <- sapply(sim.stats, function(x) x$multlik)
	 results$multinom.tailp <- length(which(results$sim.multinoms < empstats$multlik)) / Nsims
	 results$multinom.sdpd <- (results$emp.multlik - mean(results$sim.multinoms)) / sd(results$sim.multinoms)
	 }
	 
	 if("delta" %in% testStats){
	 results$emp.delta <- empstats$delta
	 results$sim.deltas <- sapply(sim.stats, function(x) x$delta)
	 results$delta.tailp <- length(which(results$sim.deltas < empstats$delta)) / Nsims
	 results$delta.sdpd <- (results$emp.delta - mean(results$sim.deltas)) /	sd(results$sim.deltas)
	 }
	 
	 if("biochemdiv" %in% testStats){
	 results$emp.biocp <- empstats$biocp
	 results$sim.biocp <- sapply(sim.stats, function(x) x$biocp)
	 results$bioc.tailp <- length(which(results$sim.biocp < empstats$biocp)) / Nsims
	 results$bioc.sdpd <- (results$emp.biocp - mean(results$sim.biocp)) / sd(results$sim.biocp)
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
	 results$meanbrsup.tailp <- length(which(results$sim.meanbrsup < empstats$brsup)) / Nsims
	 results$meanbrsup.sdpd <- (results$emp.brsup - mean(results$sim.meanbrsup)) / sd(results$sim.meanbrsup)
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
	 mahavector <- mahalanobis(all.stats.mat, colMeans(all.stats.mat), cov(all.stats.mat))
	 results$sim.maha <- mahavector[1:Nsims]
	 results$emp.maha.sdpd <- tail(mahavector, 1)
	 results$maha.tailp <- length(which(mahavector[1:Nsims] > results$emp.maha.sdpd)) / Nsims
	 }
	 
	 results$empirical.tree <- empstats$outputTree

	 if(model == "GTR+G"){
                 results$gtrMatrix <- empstats$gtrMatrix
                 results$piParams <- empstats$piParams
                 results$alphaParam <- empstats$alphaParam
         }

	 return(results)

}