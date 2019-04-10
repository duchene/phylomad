test.phylosignal <- function(sdata, format = "phylip", testType = c("locus", "genome"), aadata = F, model = "GTR+G", iqtreePath, Nsims = 100, para = F, ncore = 1, testStats = c("dnet", "dtree", "entrop", "icert", "brsup", "binp"), returnEstPhylo = F, returnAllDat = F){

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
		  write.dna(data, file = "empirical.phy")
	 } else {
	   	  ## Run ASTRAL
		  astralres <- 
		  emptre <- astralres$tree
		  write.tree(data, file = "empirical.phy")
	 }
	 if(returnEstPhylo) write.tree(emptre, file = "empirical.tre")
	 
	 # Use IQtree to calculate bipartition probabilities per branch THIS CURRENTLY ASSUMES THAT --scf and --gcf produce identical results
	 
	 emptyemptre <- emptre
	 emptyemptre$edge.length <- NULL
	 write.tree(emptyemptree, file = "empirical.empty.tre")
	 if(testType == "locus") conctype <- "scf" else conctype <- "gcf"
	 system(paste0(iqtreePath, " -t empirical.empty.tre -s empirical.phy --", conctype, " 100 --prefix emp.conc"))
	 conctab <- read.table('emp.conc.cf.stat', header=TRUE, sep = "\t")[,-7]
	 colnames(conctab)[6] <- "emp.brsup"
	 concidtr <- read.tree("emp.conc.cf.branch")
	 
	 system("rm empirical.empty.phy emp.conc.cf.tree emp.conc.log emp.conc.cf.branch emp.conc.cf.stat")
	 
	 results <- list()
	 
	 ## Calculate empirical test statistics
	 
	 if("dnet" %in% testStats) conctab$emp.dnet <- apply(conctab[,2:4]/100, 1, get.dist2net[3])
	 if("dtree" %in% testStats) conctab$emp.dtree <- apply(conctab[,2:4]/100, 1, get.dist2tr[3])
	 if("entrop" %in% testStats) conctab$emp.entrop <- apply(conctab[,2:4]/100, 1, get.quartet.entropy)
	 if("icert" %in% testStats) conctab$emp.icert <- apply(conctab[,2:4]/100, 1, get.internode.cert)
	 if("binp" %in% testStats) conctab$emp.binp <- apply(conctab[,2:5]/100, 1, get.binom.p)

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
               system(paste0(iqtreePath, " -t empirical.empty.tre -s sim.alignment.", i, ".phy --", conctype, " 100 --prefix sim.conc.", i))
               sim[[i]] <- read.table(paste0("sim.conc.", i, ".cf.stat"), header=TRUE, sep = "\t")[,-7]
               colnames(sim[[i]])[2:4] <- paste0(colnames(sim[[i]])[2:4], ".sim.", i)
	       colnames(sim[[i]])[6] <- paste0("sim.", i, ".brsup")
	       system(paste0("rm sim.alignment.", i".phy sim.conc.", i, ".cf.tree sim.conc.", i, ".log sim.conc.", i, ".cf.branch sim.conc.", i, ".cf.stat"))
	       
	 } else if(testType == "genome"){
			   




	 	   }
     
     ## Calculate test statistics for simulations. CONSIDER REMOVING THE sim LIST ALTOGETHER
     
     if("dnet" %in% testStats) sim[[i]][,paste0("sim.", i, ".dnet")] <- apply(sim[[i]][,2:4]/100, 1, get.dist2net[3])
     if("dtree" %in% testStats) sim[[i]][,paste0("sim.", i, ".dtree")] <- apply(sim[[i]][,2:4]/100, 1, get.dist2tr[3])
     if("entrop" %in% testStats) sim[[i]][,paste0("sim.", i, ".entrop")] <- apply(sim[[i]][,2:4]/100, 1, get.quartet.entropy)
     if("icert" %in% testStats) sim[[i]][,paste0("sim.", i, ".icert")] <- apply(sim[[i]][,2:4]/100, 1, get.internode.cert)
     if("binp" %in% testStats) sim[[i]][,paste0("sim.", i, ".binp")] <- apply(sim[[i]][,2:5]/100, 1, get.binom.p)    
     
     conctab <- cbind(conctab, sim[[i]])
     
     }
	 
	 ## Get P-values for test statistics
	 
	 
	 
}