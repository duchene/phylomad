source("otherScripts/R/run.gene.R")
source("otherScripts/R/get.test.statistics.R")
source("otherScripts/R/runPhyML.R")
source("otherScripts/R/clean.gene.R")
source("otherScripts/R/get.model.R")
source("testStatistics/get.chisqstat.R")
source("testStatistics/get.biodivstat.R")
source("otherScripts/R/print.bias.risk.R")

initial.dir <- getwd()

print("Functions required were loaded successfully")

if(input$treesFormat == "newick"){ trees <- read.tree(as.character(input$treesPath[1, 4])) } else if(input$treesFormat == "nexus"){ trees <- read.nexus(as.character(input$treesPath[1, 4])) } else { trees <- NULL }

if(input$treesFormat != "none" & class(trees) == "phylo"){ 
	trees <- list(trees)
	print("One tree was found as input")
}

if(input$treesFormat != "none" & length(trees) < nrow(input$dataPath)){
	trees <- rep(trees[1], nrow(input$dataPath))
	print("Since there are less trees than loci, the first tree was used for all locus analysis")
}

if(input$treesFormat != "none") class(trees) <- "multiPhylo"

selectedStats <- unlist(input$testStats)

outs <- list()

model <- paste0(input$model, input$RASmodel)

aadata <- model == "LG+G" | model == "WAG+G" | model == "JTT+G" | model == "Dayhoff+G" | model == "LG" | model == "WAG" | model == "JTT" | model == "Dayhoff"

for(j in 1:nrow(input$dataPath)){

if(input$model == "autoModel") model <- get.model(as.character(input$dataPath[j, 4]))

print("Model to be assessed was identified")

genebin <- clean.gene(sdata = as.character(input$dataPath[j, 4]), format = input$dataFormat, aadata = aadata)

print("Locus was cleaned successfully")

if(input$Ncores > 1) parallelise <- T else parallelise <- F

setwd(input$outputFolder)

print("Output folder was identified successfully")

if(input$overwirte == F && file.exists(paste0(as.character(input$dataPath[j, 1]), ".phylomad"))){
       stop("Exisitng file will not be overwritten.")
} else {
       system(paste0("mkdir ", as.character(input$dataPath[j, 1]), ".phylomad"))
       setwd(paste0(as.character(input$dataPath[j, 1]), ".phylomad"))
}

geneResults <- run.gene(sdata = genebin, format = "bin", aadata = aadata, model = model, phymlPath = phymlPath, Nsims = input$Nsims, para = parallelise, ncore = input$Ncores, testStats = selectedStats, tree = trees[j], returnSimPhylo = T, returnSimDat = T)
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

	if(any(c("chisq", "multlik", "biochemdiv", "consind", "maha") %in% selectedStats)) outs[[j]] <- cbind(outs[[j]], print.bias.risk(selectedStats, geneResults, ncol(as.matrix(genebin))))

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
	        if(aadata){
			for(i in 1:input$Nsims) write.nexus.data(geneResults$simDat[[i]], format = "protein", file = paste0("predictive.data.", i, ".nex"))
		} else {
			for(i in 1:input$Nsims) write.nexus.data(geneResults$simDat[[i]], file = paste0("predictive.data.", i, ".nex"))
		}
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
	if("chisq" %in% selectedStats) statlabels <- c(statlabels, "Chi-squared statistic")
	if("multlik" %in% selectedStats) statlabels <- c(statlabels, "Multinomial likelihood")
	if("delta" %in% selectedStats) statlabels <- c(statlabels, "Delta statistic")
	if("biochemdiv" %in% selectedStats) statlabels <- c(statlabels, "Biochemical diversity")
	if("consind" %in% selectedStats) statlabels <- c(statlabels, "Consistency Index")
	if("brsup" %in% selectedStats) statlabels <- c(statlabels, "Branch support")
	if("CIbrsup" %in% selectedStats) statlabels <- c(statlabels, "Branch support 95% interval")
	if("trlen" %in% selectedStats) statlabels <- c(statlabels, "Tree length")
	if("maha" %in% selectedStats) statlabels <- c(statlabels, "Squared Mahalanobis distance")
	if(length(empstats) == 0){
		print("Test plots cannot be returned because no test statistics were calculated.")
	} else {
	        pdf("adequacy.tests.plots.pdf", useDingbats = F, height = 5)
		for(i in 1:length(empstats)){
		      sdstat <- sd(statsmat[2:nrow(statsmat), i])
		      hist(statsmat[2:nrow(statsmat), i], xlim = c(min(statsmat[, i], na.rm = T) - sdstat, max(statsmat[, i], na.rm = T) + sdstat), xlab = statlabels[i], ylab = "Frequency of predictive simulations", main = "")
		      abline(v = empstats[i], col = "red", lwd = 3)
		}
		dev.off()
	}
}

setwd("..")

if("simple" %in% unlist(input$whatToOutput)) system(paste0("rm -r ", as.character(input$dataPath[j, 1]), ".phylomad"))

}

if(nrow(input$dataPath) > 1 && "pvals" %in% unlist(input$whatToOutput) || "simple" %in% unlist(input$whatToOutput)){
	allOutput <- do.call("rbind", outs)
	write.csv(allOutput, file = "output.all.loci.PhyloMAd.csv")
}

setwd(initial.dir)

print("Assessment completed successfully.")