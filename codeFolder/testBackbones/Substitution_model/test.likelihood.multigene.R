source("otherScripts/R/run.gene.R")
source("otherScripts/R/get.test.statistics.R")
source("otherScripts/R/runPhyML.R")
source("otherScripts/R/clean.gene.R")
source("otherScripts/R/get.model.R")
source("testStatistics/get.chisqstat.R")
source("testStatistics/get.biodivstat.R")

initial.dir <- getwd()

print("Functions required have been loaded")

if(input$treesFormat == "newick"){ trees <- read.tree(as.character(input$treesPath[1, 4])) } else if(input$treesFormat == "nexus"){ trees <- read.nexus(as.character(input$treesPath[1, 4])) } else { trees <- NULL }

if(input$treesFormat != "none" & class(trees) == "phylo"){ 
	trees <- list(trees)
	print("One tree was found as input")
}
if(input$treesFormat != "none" & length(trees) < nrow(input$dataPath)){
	trees <- rep(trees[1], nrow(input$dataPath))
	print("Since there are less trees than loci, the first tree will be used for all locus analysis")
}

if(input$treesFormat != "none") class(trees) <- "multiPhylo"

selectedStats <- unlist(input$testStats)

if("chisq" %in% selectedStats) chisqthresholds <- data.frame(seqlen = c(500, 1000, 5000), minD = c(40, 70, 367), maxD = c(135, 258, 1273))

outs <- list()

model <- paste0(input$model, input$RASmodel)

aadata <- model == "LG+G" | model == "WAG+G" | model == "JTT+G" | model == "Dayhoff+G" | model == "LG" | model == "WAG" | model == "JTT" | model == "Dayhoff"

for(j in 1:nrow(input$dataPath)){

if(input$model == "autoModel") model <- get.model(as.character(input$dataPath[j, 4]))

print("Model has been identified")

genebin <- clean.gene(sdata = as.character(input$dataPath[j, 4]), format = input$dataFormat, aadata = aadata)

print("Gene has been cleaned")

if(input$Ncores > 1) parallelise <- T else parallelise <- F

setwd(input$outputFolder)

print("Output folder has been identified")

system(paste0("mkdir ", as.character(input$dataPath[j, 1]), ".phylomad"))
setwd(paste0(as.character(input$dataPath[j, 1]), ".phylomad"))

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
	
	if("chisq" %in% selectedStats){
	       thresholdMidRisk <- predict(lm(minD ~ seqlen, data = chisqthresholds), data.frame(seqlen = ncol(as.matrix(genebin))))
	       thresholdHighRisk <- predict(lm(maxD ~ seqlen, data = chisqthresholds), data.frame(seqlen = ncol(as.matrix(genebin))))
	       if(geneResults$chisq.sdpd > thresholdHighRisk){
	       		writeLines(c("This locus is at high risk of leading to biased inferences due to comopositional heterogeneity. It is highly unadvisable to use this locus with homogeneous substitution models such as those of the GTR family.", "This advice is based on simulations in the following study:", "Duchêne, D.A., Duchêne, S., & Ho, S.Y.W. (2017). New Statistical Criteria Detect Phylogenetic Bias Caused by Compositional Heterogeneity. Molecular Biology and Evolution, 34(6), 1529-1534."), con = "locus.at.high.risk.txt")
			outs[[j]] <- cbind(outs[[j]], "high")
	       } else if(geneResults$chisq.sdpd > thresholdMidRisk){
	       	        writeLines(c("This locus is at medium risk of leading to biased inferences due to comopositional heterogeneity. These data might not lead to biased inferences when using homogeneous substitution models such as those of the GTR family. However, these inferences should be made with caution and the results of this assessment should be reported. Alternatively, a heterogeneous substitution model can be used.", "This advice is based on the following simulations study:", "Duchêne, D.A., Duchêne, S., & Ho, S.Y.W. (2017). New Statistical Criteria Detect Phylogenetic Bias Caused by Compositional Heterogeneity. Molecular Biology and Evolution, 34(6), 1529-1534."), con = "locus.at.medium.risk.txt")
			outs[[j]] <- cbind(outs[[j]], "mid")
	       } else {
			writeLines(c("This locus is at low risk of leading to biased inferences due to comopositional heterogeneity. Homogeneous substitution models such as those of the GTR family might provide reasonable results. It is also advisable to verify that other test statistics do not have extreme distances from the predictive distribution.", "This advice is based on the following simulations study:", "Duchêne, D.A., Duchêne, S., & Ho, S.Y.W. (2017). New Statistical Criteria Detect Phylogenetic Bias Caused by Compositional Heterogeneity. Molecular Biology and Evolution, 34(6), 1529-1534."), con = "locus.at.low.risk.txt")
			outs[[j]] <- cbind(outs[[j]], "low")
	       }
	       colnames(outs[[j]])[length(colnames(outs[[j]]))] <- "heterogeneity.bias.risk"
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