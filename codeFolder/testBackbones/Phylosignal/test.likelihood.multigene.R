initial.dir <- getwd()

source("otherScripts/R/test.phylosignal.R")
source("otherScripts/R/runIQtree.R")
source("otherScripts/R/clean.gene.R")
source("otherScripts/R/get.model.R")
source("testStatistics/get.phylosignal.metrics.R")
setwd("otherScripts/phybase/")
for(i in dir()) source(i)
setwd(initial.dir)

print("Functions required were loaded successfully")

selectedStats <- unlist(input$testStats)

outs <- list()

locilengths <- vector()

for(j in 1:nrow(input$dataPath)){

# FIND WAY TO DETERMINE THE TYPE OF DATA (DNA OR AA), AND SELECT MODEL WHEN DATA ARE AA

if(input$testType == "locus") model <- get.model(as.character(input$dataPath[j, 4]), format = input$dataFormat)

print("Model to be assessed was identified")

if(input$testType == "locus") analysisdata <- clean.gene(sdata = as.character(input$dataPath[j, 4]), format = input$dataFormat, aadata = F, clean = input$cleanOrNot) else analysisdata <- as.character(input$dataPath[j, 4])

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

geneResults <- test.phylosignal(sdata = analysisdata, format = if(input$testType == "locus") "bin" else input$dataFormat, testType = input$testType, aadata = F, model = model, iqtreePath = iqtreePath, astralPath = astralPath, Nsims = input$Nsims, para = parallelise, ncore = input$Ncores, testStats = selectedStats, returnEstPhylo = "phyloempres" %in% whatToOutput, returnSimulations = "simdat" %in% whatToOutput)

if("testPlots" %in% whatToOutput){
	histplotdat <- colMeans(geneResults[[1]], na.rm = T)
	names(histplotdat) <- colnames(geneResults[[1]])
	histplotdat <- histplotdat[-grep("ID|sDF|sN|gDF|gN|value|sdpd", names(histplotdat))]
	statlabels <- vector()
	if("dnet" %in% selectedStats) statlabels <- c(statlabels, "Distance to network")
	if("dtree" %in% selectedStats) statlabels <- c(statlabels, "Distance to tree")
	if("entrop" %in% selectedStats) statlabels <- c(statlabels, "Entropy")
	if("icert" %in% selectedStats) statlabels <- c(statlabels, "Internode certainty")
	if("CF" %in% selectedStats) statlabels <- c(statlabels, "Concordance factor")
	if("binp" %in% selectedStats) statlabels <- c(statlabels, "Binomial P of concordance factor")
	if("dstat" %in% selectedStats) statlabels <- c(statlabels, "D-statistic")
	if("kcstat" %in% selectedStats) statlabels <- c(statlabels, "Chifman-Kubatko statistic")
	
	if(length(selectedStats) == 0){
		print("Test plots cannot be returned because no test statistics were calculated.")
	} else {
	        pdf("tests.histograms.pdf", useDingbats = F, height = 5)
		for(i in 1:length(selectedStats)){
		      statdat <- histplotdat[grep(selectedStats[i], names(histplotdat))]
		      sdstat <- sd(statdat[2:length(statdat)])
		      hist(statdat[2:length(statdat)], xlim = c(min(statdat, na.rm = T) - sdstat, max(statdat, na.rm = T) + sdstat), xlab = statlabels[i], ylab = "Frequency of simulations", main = "")
		      abline(v = statdat[1], col = "red", lwd = 3)
		}
		dev.off()
		
		pdf("tests.summary.tree.pdf", useDingbats = F, height = if(length(geneResults[[2]]$edge.length) < 50) 15 else 30, width = if(length(geneResults[[2]]$edge.length) < 50) 30 else 60)
		par(mfrow = c(1,2))
		for(i in 1:length(selectedStats)){
		      
		      tr <- geneResults[[2]]
		      tr$edge.length <- rep(1, length(geneResults[[2]]$edge.length))
		      brpvalue <- round(geneResults[[1]][, paste0(selectedStats[i], ".p.value")], 1)
		      brsdpd <- round(geneResults[[1]][, paste0(selectedStats[i], ".sdpd")], 1)
		      names(brpvalue) <- names(brsdpd) <- geneResults[[1]][,1]
		      
		      brpvalue <- brpvalue[as.character(geneResults[[3]]$edge.length)]
		      brsdpd <- brsdpd[as.character(geneResults[[3]]$edge.length)]
		      brpvalue[is.na(brpvalue)] <- mean(brpvalue, na.rm = T)
		      brsdpd[is.na(brsdpd)] <- mean(brsdpd, na.rm = T)

		      plotBranchbyTrait(geneResults[[2]], brpvalue, mode = "edges", palette = "heat.colors", type = "fan", legend = F)
		      plotBranchbyTrait(tr, brpvalue, mode = "edges", palette = "heat.colors", type = "fan", title = paste0(statlabels[i], "\nP-value\n"))
		      plotBranchbyTrait(geneResults[[2]], brsdpd, mode = "edges", palette = "heat.colors", type = "fan", legend = F)
                      plotBranchbyTrait(tr, brsdpd, mode = "edges", palette = "heat.colors", type = "fan", title = paste0(statlabels[i], "\nz-score\n"))
		}
		dev.off()
		
		
	}
	
	if("allqp" %in% whatToOutput) write.csv(geneResults[[1]], file = "full.results.csv")

if("pvals" %in% whatToOutput){
	rownames(geneResults[[1]]) <- geneResults[[1]][,1]
	restab <- geneResults[[1]][,-grep("ID|sDF|sN|gDF|gN", colnames(geneResults[[1]]))]
	restab <- rbind(restab, colMeans(restab))
	rownames(restab)[nrow(restab)] <- "mean"
	write.csv(restab, file = "results.summary.csv")
}
}

setwd("..")

}

if(nrow(input$dataPath) > 1 & "pvals" %in% whatToOutput){
	allOutput <- do.call("rbind", outs)
	write.csv(allOutput, file = "output.all.loci.PhyloMAd.csv")
}

setwd(initial.dir)

print("Assessment completed successfully.")