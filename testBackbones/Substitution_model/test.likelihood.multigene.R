source("otherScripts/R/run.gene.R")
source("otherScripts/R/get.test.statistics.R")
source("otherScripts/R/runPhyML.R")
source("otherScripts/R/clean.gene.R")
source("otherScripts/R/get.model.R")
if("chisq" %in% unlist(input$testStats)) source("testStatistics/get.chisqstat.R")
if("biochemdiv" %in% unlist(input$testStats)) source("testStatistics/get.biodivstat.R")

print("Functions required have been loaded")

if(input$model == "autoModel") modeltested <- get.model(as.character(input$dataPath[1, 4])) else modeltested <- input$model

if(input$treesFormat == "newick"){ trees <- read.tree(as.character(input$treesPath[1, 4])) } else if(input$treesFormat == "nexus"){ trees <- read.nexus(as.character(input$treesPath[1, 4])) }

print("Model has been identified")

geneDNAbin <- clean.gene(sdata = as.character(input$dataPath[1, 4]), format = input$dataFormat)

print("Gene has been cleaned")

if(input$Ncores > 1) parallelise <- T else parallelise <- F

initial.dir <- getwd()

setwd(input$outputFolder)

print("Output folder has been identified")

if(input$treesFormat == "none"){
	geneResults <- run.gene(sdata = geneDNAbin, format = "DNAbin", model = modeltested, phymlPath = phymlPath, Nsims = input$Nsims, para = parallelise, ncore = input$Ncores, testStats = unlist(input$testStats))
} else {
	geneResults <- run.gene(sdata = geneDNAbin, format = "DNAbin", model = modeltested, phymlPath = phymlPath, Nsims = input$Nsims, para = parallelise, ncore = input$Ncores, testStats = unlist(input$testStats), tree = trees)
}

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
        system("mkdir predictive.alignments")
	setwd("predictive.alignments")
	if(input$outputFormat == "phylip"){
		for(i in 1:input$Nsims) write.dna(geneResults$simData[[i]], file = paste0("predictive.data.", i, ".phy"))
	} else if(input$outputFormat == "fasta"){
	        for(i in 1:input$Nsims) write.dna(geneResults$simData[[i]], file = paste0("predictive.data.", i, ".fasta"), format = "fasta")
	} else if(input$outputFormat == "nexus"){
	        for(i in 1:input$Nsims) write.nexus.data(geneResults$simData[[i]], file = paste0("predictive.data.", i, ".nex"))
	}
	setwd("..")
}

if("phylosimres" %in% unlist(input$whatToOutput)){
	write.tree(geneResults$simPhylos, file = "estimate.predictive.data.tre")
}

if("testPlots" %in% unlist(input$whatToOutput)){
	empstats <- unlist(geneResults[grep("emp[.]", names(geneResults))])
        simstats <- do.call(cbind, geneResults[grep("sim[.]", names(geneResults))])
	statsmat <- rbind(empstats, simstats)
	if(length(empstats) == 0){
		print("Test plots cannot be returned because no test statistics were calculated.")
	} else {
	        pdf("adequacy.tests.plots.pdf", useDingbats = F, height = 4)
		for(i in 1:length(empstats)){
		      hist(statsmat[2:nrow(statsmat), i], xlim = c(min(statsmat[, i]), max(statsmat[, i])), xlab = unlist(input$testStats)[i], main = "")
		      abline(v = empstats[i], col = "red", lwd = 3)
		}
		dev.off()
	}
}

setwd(initial.dir)