test.phylosignal <- function(sdata, format = "phylip", testType = c("locus", "genome"), aadata = F, model = "GTR+G", iqtreePath, astralPath, Nsims = 100, para = F, ncore = 1, testStats = c("dnet", "dtree", "entrop", "icert", "CF", "binp", "dstat", "kcstat"), returnEstPhylo = F, returnSimulations = F){

         if(format == "phylip"){
                  if(aadata) data <- read.aa(sdata) else data <- read.dna(sdata)
         } else if(format == "fasta"){
                  if(aadata) data <- read.aa(sdata, format = "fasta") else data <- read.dna(sdata, format = "fasta")
         } else if(format == "bin"){
                  data <- sdata
         } else if(format == "nexus" & testType == "locus"){
                  if(aadata) data <- as.AAbin(read.nexus.data(sdata)) else data <- as.DNAbin(read.nexus.data(sdata))
         } else if(format == "nexus" & testType == "genome"){
	   	  data <- read.nexus(sdata)
	 } else if(format == "newick"){
	   	  data <- read.tree(sdata)
	 }
	 
	 ## Estimate empirical locus tree or species tree
	 
	 if(testType == "locus"){
	 	  iqtreeres.emp <- runIQtree(data, format = 'bin', aadata = aadata, temp_name = "empirical", iqtreePath = iqtreePath, model = model)
		  emptre <- iqtreeres.emp$tree
		  data <- data[emptre$tip.label,]
		  write.dna(data, file = "empirical.phy")
	 } else {
	   	  ## Run ASTRAL
		  system(paste0("java -jar ", astralPath, " -i ", sdata, " -o species.tre"))
		  emptre <- read.tree("species.tre")
		  write.tree(data, file = "empirical.phy")
	 }
	 
	 # Use IQtree to calculate bipartition probabilities per branch THIS CURRENTLY ASSUMES THAT --scf and --gcf produce identical results
	 
	 emptyemptre <- emptre
	 write.tree(emptyemptre, file = "empirical.empty.tre")
	 if(testType == "locus"){
	 	  system(paste0(iqtreePath, " -t empirical.empty.tre -s empirical.phy --scf 100 --prefix emp.conc"))
	 } else {
	   	  system(paste0(iqtreePath, " -t empirical.empty.tre --gcf empirical.phy --prefix emp.conc"))
		  emptre$edge.length[is.na(emptre$edge.length)] <- 0
	 }
	 conctab <- read.table("emp.conc.cf.stat", header=TRUE, sep = "\t")
	 conctab[,7] <- round(conctab[,7], 5)
	 concidtr <- readLines("emp.conc.cf.branch")
	 concidtr <- read.tree(text = gsub(")", "):", concidtr))
	 
	 system("rm empirical empirical.phy emp.conc.cf.tree emp.conc.log emp.conc.cf.branch emp.conc.cf.stat")
	 
	 results <- list()
	 
	 ## Calculate empirical test statistics
	 
	 if("dnet" %in% testStats) conctab$emp.dnet <- apply(conctab[, 2:4]/100, 1, function(x) get.dist2net(x)[3])
	 if("dtree" %in% testStats) conctab$emp.dtree <- apply(conctab[, 2:4]/100, 1, function(x) get.dist2tr(x)[3])
	 if("entrop" %in% testStats) conctab$emp.entrop <- apply(conctab[, 2:4]/100, 1, get.quartet.entropy)
	 if("icert" %in% testStats) conctab$emp.icert <- apply(conctab[, 2:4]/100, 1, get.internode.cert)
	 if("binp" %in% testStats) conctab$emp.binp <- apply(cbind(conctab[, 2:4]/100, conctab[, 5]), 1, get.binom.p)
	 if("dstat" %in% testStats) conctab$emp.dstat <- apply(conctab[, 2:4]/100, 1, get.dstat)
	 if("kcstat" %in% testStats) conctab$emp.kcstat <- apply(conctab[, 2:4]/100, 1, get.kcstat)

	 ## Simulate data sets (alignments or genomes)
	 
     	 sim <- list()
	 for(i in 1:Nsims){
	 	   if(testType == "locus"){
	       	   l <- ncol(data)
               	   if(model == "GTR+G"){
                      rates = phangorn:::discrete.gamma(iqtreeres.emp$alphaParam, k = 4)
                      rates <- rates + 0.0001
                      sim_dat_all <- lapply(rates, function(r) simSeq(iqtreeres.emp$tree, l = round(l/4, 0), Q = iqtreeres.emp$gtrMatrix, bf = iqtreeres.emp$piParams, rate = r))
                      sim[[i]] <- as.DNAbin(c(sim_dat_all[[1]], sim_dat_all[[2]], sim_dat_all[[3]], sim_dat_all[[4]]))
               	   } else if(model == "GTR"){
                      sim[[i]] <- as.DNAbin(simSeq(iqtreeres.emp$tree, l = l, Q = iqtreeres.emp$gtrMatrix, bf = iqtreeres.emp$piParams))
               	   } else if(model == "HKY+G"){
                      rates = phangorn:::discrete.gamma(iqtreeres.emp$alphaParam, k = 4)
                      rates <- rates + 0.0001
                      sim_dat_all <- lapply(rates, function(r) simSeq(iqtreeres.emp$tree, l = round(l/4, 0), Q = c(1, 2*iqtreeres.emp$trtvRatio, 1, 1, 2*iqtreeres.emp$trtvRatio, 1), bf = iqtreeres.emp$piParams, rate = r))
                      sim[[i]] <- as.DNAbin(c(sim_dat_all[[1]], sim_dat_all[[2]], sim_dat_all[[3]], sim_dat_all[[4]]))
               	   } else if(model == "HKY"){
                      sim[[i]] <- as.DNAbin(simSeq(iqtreeres.emp$tree, l = l, Q = c(1, 2*iqtreeres.emp$trtvRatio, 1, 1, 2*iqtreeres.emp$trtvRatio, 1), bf = iqtreeres.emp$piParams))
               	   } else if(model == "JC+G"){
                      rates = phangorn:::discrete.gamma(iqtreeres.emp$alphaParam, k = 4)
                      rates <- rates + 0.0001
                      sim_dat_all <- lapply(rates, function(r) simSeq(iqtreeres.emp$tree, l = round(l/4, 0), rate = r))
                      sim[[i]] <- as.DNAbin(c(sim_dat_all[[1]], sim_dat_all[[2]], sim_dat_all[[3]], sim_dat_all[[4]]))
               	   } else if(model == "JC"){
                      sim[[i]] <- as.DNAbin(simSeq(iqtreeres.emp$tree, l = l))
               	   } else if(model == "LG+G" | model == "WAG+G" | model == "JTT+G" | model == "Dayhoff+G"){
                      simpmod <- gsub("[+]G", "", model)
                      rates = phangorn:::discrete.gamma(iqtreeres.emp$alphaParam, k = 4)
                      rates <- rates + 0.0001
                      sim_dat_all <- lapply(rates, function(r) simSeq(iqtreeres.emp$tree, type = "AA", model = simpmod, l = round(l/4, 0), rate = r))
                      sim[[i]] <- as.AAbin(c(sim_dat_all[[1]], sim_dat_all[[2]], sim_dat_all[[3]], sim_dat_all[[4]]))
               	   } else if(model == "LG" | model == "WAG" | model == "JTT" | model == "Dayhoff"){
                      sim[[i]] <- as.AAbin(simSeq(iqtreeres.emp$tree, l = l, type = "AA", model = model))
               	   }
               
		   write.dna(sim[[i]], file = paste0("sim.alignment.", i, ".phy"))
               	   system(paste0(iqtreePath, " -t empirical.empty.tre -s sim.alignment.", i, ".phy --scf 100 --prefix sim.conc.", i))
               	   sim[[i]] <- read.table(paste0("sim.conc.", i, ".cf.stat"), header=TRUE, sep = "\t")[, -c(1, 6, 7)]
               	   colnames(sim[[i]]) <- paste0(colnames(sim[[i]]), ".sim.", i)
	       	   system(paste0("rm sim.conc.", i, ".cf.tree sim.conc.", i, ".log sim.conc.", i, ".cf.branch sim.conc.", i, ".cf.stat"))
		   if(!returnSimulations) system(paste0("rm sim.alignment.", i, ".phy"))
	 } else if(testType == "genome"){
	   	   l <- length(data)
		   sim[[i]] <- list()
	 	   for(j in 1:l) sim[[i]][[j]] <- sim.coaltree.phylo(emptre)
	 	   class(sim[[i]]) <- "multiPhylo"
		   write.tree(sim[[i]], file = paste0("sim.genetrees.", i, ".phy"))
		   system(paste0(iqtreePath, " -t empirical.empty.tre --gcf sim.genetrees.", i, ".phy --prefix sim.conc.", i))
		   sim[[i]] <- read.table(paste0("sim.conc.", i, ".cf.stat"), header=TRUE, sep = "\t")[,-c(1, 6, 7)]
		   colnames(sim[[i]]) <- paste0(colnames(sim[[i]]), ".sim.", i)
               	   system(paste0("rm sim.conc.", i, ".cf.tree sim.conc.", i, ".log sim.conc.", i, ".cf.branch sim.conc.", i, ".cf.stat"))
		   if(!returnSimulations) system(paste0("rm sim.genetrees.", i, ".phy"))
	 }
     
	 ## Calculate test statistics for simulations. CONSIDER REMOVING THE sim LIST ALTOGETHER
     
	 if("dnet" %in% testStats) sim[[i]][,paste0("sim.", i, ".dnet")] <- apply(sim[[i]][,1:3]/100, 1, function(x) get.dist2net(x)[3])
     	 if("dtree" %in% testStats) sim[[i]][,paste0("sim.", i, ".dtree")] <- apply(sim[[i]][,1:3]/100, 1, function(x) get.dist2tr(x)[3])
     	 if("entrop" %in% testStats) sim[[i]][,paste0("sim.", i, ".entrop")] <- apply(sim[[i]][,1:3]/100, 1, get.quartet.entropy)
     	 if("icert" %in% testStats) sim[[i]][,paste0("sim.", i, ".icert")] <- apply(sim[[i]][,1:3]/100, 1, get.internode.cert)
     	 if("binp" %in% testStats) sim[[i]][,paste0("sim.", i, ".binp")] <- apply(cbind(sim[[i]][,1:3]/100, sim[[i]][,4]), 1, get.binom.p)
     	 if("dstat" %in% testStats) sim[[i]][,paste0("sim.", i, ".dstat")] <- apply(sim[[i]][,1:3]/100, 1, get.dstat)
     	 if("kcstat" %in% testStats) sim[[i]][,paste0("sim.", i, ".kcstat")] <- apply(sim[[i]][,1:3]/100, 1, get.kcstat)
     
	 conctab <- cbind(conctab, sim[[i]])
     
	 }
	 
	 ## Get P-values for test statistics
	
	 for(i in 1:length(testStats)){
	      stattab <- conctab[,grep(testStats[i], colnames(conctab))]
	      if(any(c("dtree", "entrop", "binp", "dstat") == testStats[i])) pvalfunc <- function(x) length(which(x[-1] > x[1]))/(length(x[-1])) else pvalfunc <- function(x) length(which(x[-1] < x[1]))/(length(x[-1]))
	      conctab[,paste0(testStats[i], ".p.value")] <- apply(stattab, 1, pvalfunc)
	      conctab[,paste0(testStats[i], ".sdpd")] <- apply(stattab, 1, function(x) (x[1] - mean(x[-1], na.rm = T)) / sd(x[-1], na.rm = T))
	 }
	
	 #if(!returnAllDat) conctab[,-grep("ID|sCF|sDF|sN|gCF|gDF|gN", colnames(conctab))]
	 results <- list(testsTable = conctab, treeEstimate = emptre, IDtree = concidtr)
	 system("rm empirical.empty.tre")
	 if(testType == "locus" & returnEstPhylo) write.tree(emptre, file = "empirical.estimated.tre")
	 return(results)
}