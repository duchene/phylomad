source("otherScripts/R/run.gene.clock.R")
source("otherScripts/R/get.test.statistics.R")
source("otherScripts/R/runPhyML.R")
source("otherScripts/R/getRatogs.R")
if("chisq" %in% unlist(input$testStats)) source("testStatistics/get.chisqstat.R")
if("biochemdiv" %in% unlist(input$testStats)) source("testStatistics/get.biodivstat.R")
if("df" %in% unlist(input$testStats)) source("testStatistics/get.df.R")
if("stemmystat" %in% unlist(input$testStats)) source("testStatistics/stemmy.R")


initial.dir <- getwd()

print("Functions required have been loaded")

trees <- read.nexus(as.character(input$treesPath[1, 4]))

for(j in 1:nrow(input$dataPath)){

if(input$Ncores > 1) parallelise <- T else parallelise <- F

setwd(input$outputFolder)

print("Output folder has been identified")

system(paste0("mkdir ", as.character(input$dataPath[j, 1]), ".phylomad.clock"))
setwd(paste0(as.character(input$dataPath[j, 1]), ".phylomad.clock"))

selectedStats <- unlist(input$testStats)

if(nrow(input$treesPath) == 1) treesPath <- as.character(input$treesPath[1, 4]) else treesPath <- as.character(input$treesPath[j, 4])
if(nrow(input$posteriorPath) == 1) postPath <- as.character(input$posteriorPath[1, 4]) else postPath <- as.character(input$posteriorPath[j, 4])

geneResults <- run.gene.clock(sdata = as.character(input$dataPath[j, 4]), format = input$dataFormat, treesFile = treesPath, logFile = postPath, burninpercentage = input$burnin, phymlPath = phymlPath, Nsims = input$Nsims, para = parallelise, ncore = input$Ncores, testStats = selectedStats, returnSimPhylo = T, returnSimDat = T)

if("pvals" %in% unlist(input$whatToOutput)){
	out <- rbind(unlist(geneResults[grep("[.]tailp", names(geneResults))]), unlist(geneResults[grep("emp[.]", names(geneResults))]), unlist(geneResults[grep("[.]sdpd", names(geneResults))]))
	if(length(out) == 0){
	       print("P-values cannot be returned because no test statistics were calculated.")
	} else {
	       colnames(out) <- unlist(input$testStats)
	       rownames(out) <- c("Tail area probability", "Empirical test statistic", "Standard deviations from simulated distribution")
	       write.csv(out, file = "output.pvals.PhyloMAd.csv")
	}
}

if("phyloempres" %in% unlist(input$whatToOutput)){
	write.tree(geneResults$empirical.tree, file = "estimate.empirical.data.tre")
}

if("simdat" %in% unlist(input$whatToOutput)){
	if(input$outputFormat == "phylip"){
		for(i in 1:input$Nsims) write.dna(geneResults$simDat[[i]], file = paste0("predictive.data.", i, ".phy"))
	} else if(input$outputFormat == "fasta"){
	        for(i in 1:input$Nsims) write.dna(geneResults$simDat[[i]], file = paste0("predictive.data.", i, ".fasta"), format = "fasta")
	} else if(input$outputFormat == "nexus"){
	        for(i in 1:input$Nsims) write.nexus.data(geneResults$simDat[[i]], file = paste0("predictive.data.", i, ".nex"))
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
	}
}

setwd("..")

}

setwd(initial.dir)