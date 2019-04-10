source("otherScripts/R/test.phylosignal.R")
source("otherScripts/R/runIQtree.R")
source("otherScripts/R/clean.gene.R")
source("otherScripts/R/get.model.R")
source("testStatistics/get.phylosignal.metrics.R")

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

locilengths <- vector()

for(j in 1:nrow(input$dataPath)){

if(input$model == "autoModel") model <- get.model(as.character(input$dataPath[j, 4]))

print("Model to be assessed was identified")

genebin <- clean.gene(sdata = as.character(input$dataPath[j, 4]), format = input$dataFormat, aadata = aadata, clean = input$cleanOrNot)

print("Locus was cleaned successfully")

if(input$Ncores > 1) parallelise <- T else parallelise <- F

setwd(input$outputFolder)

print("Output folder was identified successfully")

if(!input$overwrite && file.exists(paste0(as.character(input$dataPath[j, 1]), ".phylomad"))){
       stop("Exisitng files will not be overwritten. Aborting.")
} else {
       system(paste0("mkdir ", as.character(input$dataPath[j, 1]), ".phylomad"))
       setwd(paste0(as.character(input$dataPath[j, 1]), ".phylomad"))
}

whatToOutput <- unlist(input$whatToOutput)

geneResults <- run.gene(sdata = genebin, format = "bin", aadata = aadata, model = model, iqtreePath = iqtreePath, Nsims = input$Nsims, para = parallelise, ncore = input$Ncores, testStats = selectedStats, tree = trees[j], returnSimPhylo = T, returnSimDat = T)

if("pvals" %in% whatToOutput || "simple" %in% whatToOutput){
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

	if(any(c("chisq", "multlik", "biochemdiv", "consind", "maha") %in% selectedStats)){
		locilengths[j] <- ncol(as.matrix(genebin))
		biasrisk <- print.bias.risk(selectedStats, geneResults, locilengths[j])
		outs[[j]] <- cbind(outs[[j]], biasrisk[[1]])
	}

}

if("phyloempres" %in% whatToOutput){
	write.tree(geneResults$empirical.tree, file = "estimate.empirical.data.tre")
}

if("simdat" %in% whatToOutput){
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

if("phylosimres" %in% whatToOutput){
	write.tree(geneResults$simPhylos, file = "estimate.predictive.data.tre")
}

if("testPlots" %in% whatToOutput){
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

if("simple" %in% whatToOutput) system(paste0("rm -r ", as.character(input$dataPath[j, 1]), ".phylomad"))

}

if(nrow(input$dataPath) > 1 && "pvals" %in% whatToOutput || "simple" %in% whatToOutput){
	allOutput <- do.call("rbind", outs)
	write.csv(allOutput, file = "output.all.loci.PhyloMAd.csv")
}

if("multiTestPlots" %in% whatToOutput){
	pdf("multi.locus.results.plots.pdf", height = 4, width = 4, useDingbats = F)
	if("chisq" %in% selectedStats && !any(allOutput[,"chisq.stdev.from.pred.dist"] == Inf) && !any(is.na(allOutput[,"chisq.stdev.from.pred.dist"]))){ 
		plot(as.numeric(allOutput[,"chisq.stdev.from.pred.dist"]), log(locilengths), main = "Chi-squared statistic", xlab = "Standard deviations from mean\nof predictive distribution (SDPD)", ylab = "Log locus length (number of sites)", pch = 19, xlim = c(if(min(as.numeric(allOutput[,"chisq.stdev.from.pred.dist"])) < 0) min(as.numeric(allOutput[,"chisq.stdev.from.pred.dist"])) else 0, max(as.numeric(allOutput[,"chisq.stdev.from.pred.dist"]))), ylim = c(0, if(max(locilengths) < 500) 500 else if(max(locilengths) < 5000) 5000 else max(locilengths)))
		abline(lm(biasrisk[[2]][,"seqlen.chisq"] ~ biasrisk[[2]][,"minD.chisq"]), lwd = 2, col = "orange")
		abline(lm(biasrisk[[2]][,"seqlen.chisq"] ~ biasrisk[[2]][,"maxD.chisq"]), lwd = 2, col = "red")
		points(biasrisk[[2]][,"minD.chisq"], biasrisk[[2]][,"seqlen.chisq"], pch = 19, col = "orange")
		points(biasrisk[[2]][,"maxD.chisq"], biasrisk[[2]][,"seqlen.chisq"], pch = 19, col = "red")
	}
	if("multlik" %in% selectedStats && !any(allOutput[,"multlik.stdev.from.pred.dist"] == Inf) && !any(is.na(allOutput[,"multlik.stdev.from.pred.dist"]))){
		plot(as.numeric(allOutput[,"multlik.stdev.from.pred.dist"]), log(locilengths), main = "Multinomial likelihood statistic", xlab = "Standard deviations from mean\nof predictive distribution (SDPD)", ylab = "Log locus length (number of sites)", pch = 19, xlim = c(if(min(as.numeric(allOutput[,"multlik.stdev.from.pred.dist"])) < 0) min(as.numeric(allOutput[,"multlik.stdev.from.pred.dist"])) else 0, max(as.numeric(allOutput[,"multlik.stdev.from.pred.dist"]))), ylim = c(0, if(max(locilengths) < 500) log(500) else if(max(locilengths) < 5000) log(5000) else log(max(locilengths))))
		abline(lm(log(biasrisk[[2]][,"seqlen"]) ~ biasrisk[[2]][,"minD.multlik"]), lwd = 2, col = "orange")
		abline(lm(log(biasrisk[[2]][,"seqlen"]) ~ biasrisk[[2]][,"maxD.multlik"]), lwd = 2, col = "red")
		points(biasrisk[[2]][,"minD.multlik"], log(biasrisk[[2]][,"seqlen"]), pch = 19, col = "orange")
		points(biasrisk[[2]][,"maxD.multlik"], log(biasrisk[[2]][,"seqlen"]), pch = 19, col = "red")
	}
	if("biochemdiv" %in% selectedStats && !any(allOutput[,"biochemdiv.stdev.from.pred.dist"] == Inf) && !any(is.na(allOutput[,"biochemdiv.stdev.from.pred.dist"]))){
		plot(as.numeric(allOutput[,"biochemdiv.stdev.from.pred.dist"]), log(locilengths), main = "Biochemical diversity statistic", xlab = "Standard deviations from mean\nof predictive distribution (SDPD)", ylab = "Log locus length (number of sites)", pch = 19, xlim = c(if(max(as.numeric(allOutput[,"biochemdiv.stdev.from.pred.dist"])) > 0) max(as.numeric(allOutput[,"biochemdiv.stdev.from.pred.dist"])) else 0, min(as.numeric(allOutput[,"biochemdiv.stdev.from.pred.dist"]))), ylim = c(0, if(max(locilengths) < 500) log(500) else if(max(locilengths) < 5000) log(5000) else log(max(locilengths))))
		abline(lm(log(biasrisk[[2]][,"seqlen"]) ~ biasrisk[[2]][,"minD.biochemdiv"]), lwd = 2, col = "orange")
		abline(lm(log(biasrisk[[2]][,"seqlen"]) ~ biasrisk[[2]][,"maxD.biochemdiv"]), lwd = 2, col = "red")
		points(biasrisk[[2]][,"minD.biochemdiv"], log(biasrisk[[2]][,"seqlen"]), pch = 19, col = "orange")
		points(biasrisk[[2]][,"maxD.biochemdiv"], log(biasrisk[[2]][,"seqlen"]), pch = 19, col = "red")
	}
	if("consind" %in% selectedStats && !any(allOutput[,"consind.stdev.from.pred.dist"] == Inf) && !any(is.na(allOutput[,"consind.stdev.from.pred.dist"]))){
		plot(as.numeric(allOutput[,"consind.stdev.from.pred.dist"]), log(locilengths), main = "Consistency Index statistic", xlab = "Standard deviations from mean\nof predictive distribution (SDPD)", ylab = "Log locus length (number of sites)", pch = 19, xlim = c(if(max(as.numeric(allOutput[,"consind.stdev.from.pred.dist"])) > 0) max(as.numeric(allOutput[,"consind.stdev.from.pred.dist"])) else 0, min(as.numeric(allOutput[,"consind.stdev.from.pred.dist"]))), ylim = c(0, if(max(locilengths) < 500) log(500) else if(max(locilengths) < 5000) log(5000) else log(max(locilengths))))
		abline(lm(log(biasrisk[[2]][,"seqlen"]) ~ biasrisk[[2]][,"minD.consind"]), lwd = 2, col = "orange")
		abline(lm(log(biasrisk[[2]][,"seqlen"]) ~ biasrisk[[2]][,"maxD.consind"]), lwd = 2, col = "red")
		points(biasrisk[[2]][,"minD.consind"], log(biasrisk[[2]][,"seqlen"]), pch = 19, col = "orange")
		points(biasrisk[[2]][,"maxD.consind"], log(biasrisk[[2]][,"seqlen"]), pch = 19, col = "red")
	}
	if("maha" %in% selectedStats && !any(allOutput[,"maha.stdev.from.pred.dist"] == Inf) && !any(is.na(allOutput[,"maha.stdev.from.pred.dist"]))){
		plot(as.numeric(allOutput[,"maha.stdev.from.pred.dist"]), log(locilengths), main = "Mahalanobis distance statistic", xlab = "Standard deviations from mean\nof predictive distribution (SDPD)", ylab = "Log locus length (number of sites)", pch = 19, xlim = c(0, max(as.numeric(allOutput[,"maha.stdev.from.pred.dist"]))), ylim = c(0, if(max(locilengths) < 500) log(500) else if(max(locilengths) < 5000) log(5000) else log(max(locilengths))))
		abline(lm(log(biasrisk[[2]][,"seqlen"]) ~ biasrisk[[2]][,"minD.maha"]), lwd = 2, col = "orange")
		points(biasrisk[[2]][,"minD.maha"], log(biasrisk[[2]][,"seqlen"]), pch = 19, col = "orange")
	}
	dev.off()
}

setwd(initial.dir)

print("Assessment completed successfully.")