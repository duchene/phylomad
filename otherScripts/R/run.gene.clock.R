run.gene.clock <- function(sdata, treesFile, logFile, burninpercentage, format = "phylip", phymlPath, Nsims = 100, para = F, ncore = 1, testStats = c("imbal", "stemmystat", "df", "trlen"), returnSimPhylo = F, returnSimDat = F){
	 
	 # Load data
	 
	 if(format == "phylip"){
                  data <- read.dna(sdata)
         } else if(format == "fasta"){
                  data <- read.dna(sdata, format = "fasta")
         } else if(format == "DNAbin"){
                  data <- sdata
         } else if(format == "nexus"){
	   	  data <- as.DNAbin(read.nexus.data(sdata))
	 }
	 
	 trees <- read.nexus(treesFile)
         logdat <- read.table(logFile, header = T, comment = "#")
         burninsamples <- 1:(round(length(trees) * (burninpercentage/100)))
	 trees <- trees[-burninsamples]
	 logdat <- logdat[-burninsamples,]
	 samp <- sample(1:length(trees), Nsims)
         trees <- trees[samp]
         logdat <- logdat[samp,]
	 
	 if("ucldMean" %in% colnames(logdat) | "meanClockRate" %in% colnames(logdat)){
                  ratogs <- getRatogs(trees.file)[samp]
         }

	 # Simulate data sets.

	 l <- ncol(as.matrix(data))
	 sim <- list()
	 for(i in 1:Nsims){
	 
               tr <- trees[[i]]
               if("clockRate" %in% colnames(logdat)){
                      tr$edge.length <- tr$edge.length * logdat[i, "clockRate"]
                      sim[[i]] <- list(phylogram = tr)
               } else if("ucldMean" %in% colnames(logdat) | "meanClockRate" %in% colnames(logdat)){
                      trr <- ratogs[[i]]
                      trr$edge.length <- trr$edge.length * tr$edge.length
                      sim[[i]] <- list(phylogram = trr)
               }
	       
	       if(all(c("rateAC", "rateAG", "rateAT", "rateCG", "rateGT") %in% colnames(logdat))){
                      # GENERAL TIME REVERSIBLE (GTR)

                      basef <- c(logdat$freqParameter.1[i], logdat$freqParameter.2[i], logdat$freqParameter.3[i], logdat$freqParameter.4[i])
                      qmat <- c(logdat$rateAC[i], logdat$rateAG[i], logdat$rateAT[i], logdat$rateCG[i], 1, logdat$rateGT[i])

                      if("gammaShape" %in% colnames(logdat)){
		           model <- "GTR+G"
                           rates = phangorn:::discrete.gamma(logdat$gammaShape[i], k = 4)
			   rates <- rates + 0.0001
                           sim_dat_all<- lapply(rates, function(r) simSeq(sim[[i]][[1]], l = round(l/4, 0), Q = qmat, bf = basef, rate = r))
                           sim[[i]]$alignment <- as.DNAbin(c(sim_dat_all[[1]], sim_dat_all[[2]], sim_dat_all[[3]], sim_dat_all[[4]]))
                      } else {
		           model <- "GTR"
                           sim[[i]]$alignment <- as.DNAbin(simSeq(sim[[i]][[1]], Q = qmat, bf = basef, l = l))
                      }

               } else if("kappa" %in% colnames(logdat)){
                      # HASEGAWA-KISHINO-YANO (HKY)

                      basef <- c(logdat$freqParameter.1[i], logdat$freqParameter.2[i], logdat$freqParameter.3[i], logdat$freqParameter.4[i])
                      qmat <- c(1, 2*logdat$kappa[i], 1, 1, 2*logdat$kappa[i], 1)

                      if("gammaShape" %in% colnames(logdat)){
                           model <- "HKY+G"
			   rates = phangorn:::discrete.gamma(logdat$gammaShape[i], k = 4)
			   rates <- rates + 0.0001
                           sim_dat_all<- lapply(rates, function(r) simSeq(sim[[i]][[1]], l = round(l/4, 0), Q = qmat, bf = basef, rate = r))
                           sim[[i]]$alignment <- as.DNAbin(c(sim_dat_all[[1]], sim_dat_all[[2]], sim_dat_all[[3]], sim_dat_all[[4]]))
                      } else {
                           model <- "HKY"
			   sim[[i]]$alignment <- as.DNAbin(simSeq(sim[[i]][[1]], Q = qmat, bf = basef, l = l))
                      }

               } else {
                      # JUKES-CANTOR (JC)
		      
		      if("gammaShape" %in% colnames(logdat)){
		           model <- "JC+G"
			   rates = phangorn:::discrete.gamma(logdat$gammaShape[i], k = 4)
                           rates <- rates + 0.0001
                           sim_dat_all<- lapply(rates, function(r) simSeq(sim[[i]][[1]], l = round(l/4, 0), rate = r))
                           sim[[i]]$alignment <- as.DNAbin(c(sim_dat_all[[1]], sim_dat_all[[2]], sim_dat_all[[3]], sim_dat_all[[4]]))
		      } else {
		      	   model <- "JC"
			   sim[[i]]$alignment <- as.DNAbin(simSeq(sim[[i]][[1]], l = l))
		      }

               }
	 
	 }
	 
	 # Get test statistics for empirical data
	 
	 
	 empstats <- get.test.statistics(data, format = "DNAbin", geneName = "empirical.alignment.phy", phymlPath = phymlPath, model = model, stats = testStats, tree = trees[[i]])
         system("rm empirical.alignment.phy")

	 # Get test statistics for simulations

	 if(!para){
	   sim.stats <- list()
	 
	   for(i in 1:Nsims){	       
	       sim.stats[[i]] <- get.test.statistics(sim[[i]][[2]], format = "DNAbin", geneName = paste0("sim.data.", i), phymlPath = phymlPath, model = model, stats = testStats, tree = trees[[i]])
	       system(paste0("rm ", paste0("sim.data.", i)))
	   }
	   
	 } else {
	   ### START PARALLEL COMPUTING
	   print('Enter parallel computing')
	   require(foreach)
	   require(doParallel)
		
	   runSim <- function(i){
	     tRep <- get.test.statistics(sim[[i]], format = "DNAbin", geneName = paste0("sim.data.", i), phymlPath = phymlPath, model = model, stats = testStats, tree = trees[[i]])
             system(paste0("rm ", paste0("sim.data.", i)))
	     return(tRep)		
	   }	  
	   cl <- makeCluster(ncore)
	   registerDoParallel(cl)
	   simReps <- foreach(x = 1:Nsims, .packages = c('phangorn', 'ape'), .export = c('get.test.statistics', 'runPhyML', 'get.chisqstat', 'get.biodivstat')) %dopar% runSim(x)
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
	 failstats <- vector()
	 if(length(all.emp.stats) == 0){
	 print("No test statistics were calculated.")
	 } else {
	 for(i in 1:ncol(all.stats.mat)){
	        if(length(unique(all.sim.stats[,i])) == 1){
			print(paste(colnames(all.stats.mat)[i], "will not be included in the mahalanobis calculation and is likely to be unreliable. This is possibly because the values for all simulations are the same."))
	 		failstats <- c(failstats, i)
	 	}
	 }
	 if(length(failstats) > 0){
	 	all.sim.stats <- all.sim.stats[,-failstats]
           	all.stats.mat <- all.stats.mat[,-failstats]
	 }
	 #print(all.stats.mat)
	 if(length(all.stats.mat) >= 2*Nsims){
	 	mahavector <- mahalanobis(all.stats.mat, colMeans(all.stats.mat), cov(all.stats.mat))
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