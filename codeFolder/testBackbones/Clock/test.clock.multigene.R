source("otherScripts/R/run.gene.clock.R")
source("otherScripts/R/get.test.statistics.R")
source("otherScripts/R/runPhyML.R")
source("otherScripts/R/getRatogs.R")
source("testStatistics/get.df.R")
source("testStatistics/stemmy.R")

initial.dir <- getwd()

print("Functions required were loaded successfully")

trees <- read.nexus(as.character(input$treesPath[1, 4]))

selectedStats <- unlist(input$testStats)

outs <- list()

for(j in 1:nrow(input$dataPath)){

if(input$Ncores > 1) parallelise <- T else parallelise <- F

setwd(input$outputFolder)

print("Output folder was identified")

if(input$overwrite == F && file.exists(paste0(as.character(input$dataPath[j, 1]), ".phylomad.clock"))){
       stop("Existing file will not be overwritten")
} else {
       system(paste0("mkdir ", as.character(input$dataPath[j, 1]), ".phylomad.clock"))
       setwd(paste0(as.character(input$dataPath[j, 1]), ".phylomad.clock"))
}

if(nrow(input$treesPath) == 1) treesPath <- as.character(input$treesPath[1, 4]) else treesPath <- as.character(input$treesPath[j, 4])
if(nrow(input$posteriorPath) == 1) postPath <- as.character(input$posteriorPath[1, 4]) else postPath <- as.character(input$posteriorPath[j, 4])

geneResults <- run.gene.clock(sdata = as.character(input$dataPath[j, 4]), format = input$dataFormat, treesFile = treesPath, logFile = postPath, burninpercentage = input$burnin, phymlPath = phymlPath, Nsims = input$Nsims, para = parallelise, ncore = input$Ncores, testStats = selectedStats, returnSimPhylo = T, returnSimDat = T)

if("pvals" %in% unlist(input$whatToOutput) || "simple" %in% unlist(input$whatToOutput)){
	geneResMat <- matrix(NA, nrow = 3, ncol = length(selectedStats))
	for(i in 1:length(selectedStats)){
	       geneResMat[1, i] <- geneResults[[paste0(selectedStats[i], ".tailp")]]
	       geneResMat[2, i] <- geneResults[[paste0("emp.", selectedStats[i])]]
	       geneResMat[3, i] <- geneResults[[paste0(selectedStats[i], ".sdpd")]]
	}
	outs[[j]] <- geneResMat
	if(length(outs[[j]]) == 0){
	       print("P-values cannot be returned because no test statistics were calculated.")
	} else {
	       colnames(outs[[j]]) <- unlist(input$testStats)
	       rownames(outs[[j]]) <- c("Tail area probability", "Empirical test statistic", "Standard deviations from simulated distribution")
	       write.csv(outs[[j]], file = "output.pvals.PhyloMAd.csv")
	}
	resvector <- matrix(as.vector(t(outs[[j]])), nrow = 1)
        rownames(resvector) <- as.character(input$dataPath[j, 1])
        colnames(resvector) <- c(paste0(colnames(outs[[j]]), ".tail.area.p"), paste0(colnames(outs[[j]]), ".empirical.statistic"), paste0(colnames(outs[[j]]), ".stdev.from.pred.dist"))
	outs[[j]] <- resvector
	
}

if("phyloempres" %in% unlist(input$whatToOutput)){
	if("aindex" %in% unlist(input$testStats)){
	       branchwise.assessment <- rbind(geneResults$aindex.allp, geneResults$aindex.sds)
	       rownames(branchwise.assessment) <- c("Branch-wise P-values", "Branch-wise SDPD")
	       write.csv(branchwise.assessment, file = "Branch-wise.assessment.results.csv")
	}
	write.tree(geneResults$empirical.tree, file = "estimate.empirical.data.clock.free.tre")
}

if("simdat" %in% unlist(input$whatToOutput)){
	if(input$outputFormat == "phylip"){
		for(i in 1:input$Nsims) write.dna(geneResults$simDat[[i]]$alignment, file = paste0("predictive.data.", i, ".phy"))
	} else if(input$outputFormat == "fasta"){
	        for(i in 1:input$Nsims) write.dna(geneResults$simDat[[i]]$alignment, file = paste0("predictive.data.", i, ".fasta"), format = "fasta")
	} else if(input$outputFormat == "nexus"){
	        for(i in 1:input$Nsims) write.nexus.data(geneResults$simDat[[i]]$alignment, file = paste0("predictive.data.", i, ".nex"))
	}
}

if("phylosimres" %in% unlist(input$whatToOutput)){
	write.tree(geneResults$simPhylos, file = "estimate.predictive.data.tre")
}

if("testPlots" %in% unlist(input$whatToOutput)){
	empstats <- unlist(geneResults[grep("emp[.]", names(geneResults))])
        simstats <- do.call(cbind, geneResults[grep("sim[.]", names(geneResults))])
	statsmat <- rbind(empstats, simstats)
	statlabels <- vector()
	if("imbal" %in% selectedStats) statlabels <- c(statlabels, "Imbalance (Colless index)")
	if("stemmystat" %in% selectedStats) statlabels <- c(statlabels, "Stemminess")
	if("df" %in% selectedStats) statlabels <- c(statlabels, "Df neutrality statistic")
	if("aindex" %in% selectedStats) statlabels <- c(statlabels, "A index")
	if("trlen" %in% selectedStats) statlabels <- c(statlabels, "Tree length")
	if("maha" %in% selectedStats) statlabels <- c(statlabels, "Squared Mahalanobis distance")
	
	if(length(empstats) == 0){
		print("Test plots cannot be returned because no test statistics were calculated.")
	} else {
	        pdf("adequacy.tests.plots.pdf", useDingbats = F, height = 5)
		for(i in 1:length(empstats)){
		      sdstat <- sd(statsmat[2:nrow(statsmat), i])
		      hist(statsmat[2:nrow(statsmat), i], xlim = c(min(statsmat[, i]) - sdstat, max(statsmat[, i]) + sdstat), xlab = statlabels[i], ylab = "Frequency of predictive simulations", main = "")
		      abline(v = empstats[i], col = "red", lwd = 3)
		}
		dev.off()
		
		if("aindex" %in% unlist(input$testStats)){
		      pdf("adequacy.branch-wise.tests.plots.pdf", useDingbats = F, height = 5 + (Ntip(geneResults$empirical.tree) * 0.2))
                      treetoprint <- geneResults$empirical.tree
                      treetoprint$edge.length <- NULL
                      brPal <- colorRampPalette(c('blue', 'green', 'yellow','red'))
                      Aindexcol <- brPal(5)[as.numeric(cut(geneResults$aindex.allp, breaks = 5))]
                      sdcol <- brPal(5)[as.numeric(cut(geneResults$aindex.sds, breaks = 5))]
                      plot(treetoprint, col = Aindexcol, main = "Branches labelled by P-value")
                      edgelabels(round(geneResults$aindex.allp, 2))
                      plot(treetoprint, col = sdcol, main = "Branches labelled by SDPD")
                      edgelabels(round(geneResults$aindex.sds, 2))
		      dev.off()
                }
	}
	
}

setwd("..")

if("simple" %in% unlist(input$whatToOutput)) system(paste0("rm -r ", as.character(input$dataPath[j, 1]), ".phylomad.clock"))

}

if(nrow(input$dataPath) > 1 && "pvals" %in% unlist(input$whatToOutput) || "simple" %in% unlist(input$whatToOutput)){
        allOutput <- do.call("rbind", outs)
        write.csv(allOutput, file = "output.all.loci.clock.PhyloMAd.csv")
}

setwd(initial.dir)