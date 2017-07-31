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
                  ratogs <- getRatogs(treesFile)[samp]
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
	 
	 
	 empstats <- get.test.statistics(data, format = "bin", geneName = "empirical.alignment.phy", phymlPath = phymlPath, model = model, stats = testStats, tree = trees[[i]])
         system("rm empirical.alignment.phy")

	 # Get test statistics for simulations

	 if(!para){
	   sim.stats <- list()
	 
	   for(i in 1:Nsims){	       
	       sim.stats[[i]] <- get.test.statistics(sim[[i]]$alignment, format = "bin", geneName = paste0("sim.data.", i), phymlPath = phymlPath, model = model, stats = testStats, tree = trees[[i]])
	       system(paste0("rm ", paste0("sim.data.", i)))
	   }
	   
	 } else {
	   ### START PARALLEL COMPUTING
	   print('Parallel computing started')
	   require(foreach)
	   require(doParallel)
		
	   runSim <- function(i){
	     tRep <- get.test.statistics(sim[[i]]$alignment, format = "bin", geneName = paste0("sim.data.", i), phymlPath = phymlPath, model = model, stats = testStats, tree = trees[[i]])
             system(paste0("rm ", paste0("sim.data.", i)))
	     return(tRep)		
	   }	  
	   cl <- makeCluster(ncore)
	   registerDoParallel(cl)
	   simReps <- foreach(x = 1:Nsims, .packages = c('phangorn', 'ape', 'apTreeshape'), .export = c('get.test.statistics', 'runPhyML', 'get.df', 'stemmy')) %dopar% runSim(x)
	   sim.stats <- simReps 
	   stopCluster(cl)
	   print("Parallel computing ended successfully")
	  }
	  ### END PARALLEL COMPUTING

	 # Get P-values for test statistics.
	 
	 results <- list()
	 
	 if("trlen" %in% testStats){
	 results$emp.trlen <- empstats$trlen
	 results$sim.trlens <- sapply(sim.stats, function(x) x$trlen)
	 results$trlen.tailp <- length(which(results$sim.trlens < empstats$trlen)) / Nsims
	 results$trlen.sdpd <- (results$emp.trlen - mean(results$sim.trlens)) /	sd(results$sim.trlens)
	 }

	 if("imbal" %in% testStats){
	 results$emp.imbal <- empstats$imbal
	 results$sim.imbals <- sapply(sim.stats, function(x) x$imbal)
	 results$imbal.tailp <- length(which(results$sim.imbals < empstats$imbal)) / Nsims
	 results$imbal.sdpd <- (results$emp.imbal - mean(results$sim.imbals)) / sd(results$sim.imbals)
	 }

	 if("stemmystat" %in% testStats){
	 results$emp.stemmystat <- empstats$stemmystat
	 results$sim.stemmystats <- sapply(sim.stats, function(x) x$stemmystat)
	 results$stemmystat.tailp <- length(which(results$sim.stemmystats < empstats$stemmystat)) / Nsims
	 results$stemmystat.sdpd <- (results$emp.stemmystat - mean(results$sim.stemmystats, na.rm = T)) / sd(results$sim.stemmystats, na.rm = T)
	 }
	 
	 if("df" %in% testStats){
	 results$emp.df <- empstats$df
	 results$sim.dfs <- sapply(sim.stats, function(x) x$df)
	 results$df.tailp <- length(which(results$sim.dfs < empstats$df)) / Nsims
	 results$df.sdpd <- (results$emp.df - mean(results$sim.dfs)) /	sd(results$sim.dfs)
	 }
	 
	 if("aindex" %in% testStats){
	 results$realdat.aindex <- empstats$aindex
	 results$emp.aindex <- median(empstats$aindex)
	 results$all.predblens <- lapply(sim.stats, function(x) x$aindex)
	 results$sim.aindexs <- sapply(results$all.predblens, median)
	 results$aindex.allp <- vector()
	 results$aindex.sds <- vector()
	 for(i in 1:length(empstats$aindex)){
	        branchpredlens <- sapply(results$all.predblens, function(x) x[i])
	 	results$aindex.allp[i] <- length(which(branchpredlens < empstats$aindex[i])) / Nsims
	 	results$aindex.sds[i] <- (results$realdat.aindex[i] - mean(branchpredlens)) / sd(branchpredlens)
	 }
	 results$aindex.tailp <- length(which(results$aindex.allp > 0.05)) / length(results$aindex.allp)
	 results$aindex.sdpd <- mean(results$aindex.sds)
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
			print(paste(colnames(all.stats.mat)[i], "was not included in the mahalanobis calculation and is likely to be unreliable. This is possibly because the values for all simulations are the same."))
	 		failstats <- c(failstats, i)
	 	}
	 }
	 if(length(failstats) > 0){
	 	all.sim.stats <- all.sim.stats[,-failstats]
           	all.stats.mat <- all.stats.mat[,-failstats]
	 }
	 
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