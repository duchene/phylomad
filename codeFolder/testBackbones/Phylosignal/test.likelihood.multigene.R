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
	empstats <- unlist(geneResults[grep("emp[.]", names(geneResults))])
        simstats <- do.call(cbind, geneResults[grep("sim[.]", names(geneResults))])
	statsmat <- rbind(empstats, simstats)
	statlabels <- vector()
	if("dnet" %in% selectedStats) statlabels <- c(statlabels, "Distance to network")
	if("dtree" %in% selectedStats) statlabels <- c(statlabels, "Distance to tree")
	if("entrop" %in% selectedStats) statlabels <- c(statlabels, "Entropy")
	if("icert" %in% selectedStats) statlabels <- c(statlabels, "Internode certainty")
	if("CF" %in% selectedStats) statlabels <- c(statlabels, "Concordance factor")
	if("binp" %in% selectedStats) statlabels <- c(statlabels, "Binomial P of concordance factor")
	if("dstat" %in% selectedStats) statlabels <- c(statlabels, "D-statistic")
	if("kcstat" %in% selectedStats) statlabels <- c(statlabels, "Chifman-Kubatko statistic")
	
	if(length(empstats) == 0){
		print("Test plots cannot be returned because no test statistics were calculated.")
	} else {
	        pdf("tests.histograms.pdf", useDingbats = F, height = 5)
		for(i in 1:length(empstats)){
		      sdstat <- sd(statsmat[2:nrow(statsmat), i])
		      hist(statsmat[2:nrow(statsmat), i], xlim = c(min(statsmat[, i], na.rm = T) - sdstat, max(statsmat[, i], na.rm = T) + sdstat), xlab = statlabels[i], ylab = "Frequency of predictive simulations", main = "")
		      abline(v = empstats[i], col = "red", lwd = 3)
		}
		dev.off()
		
		pdf("tests.summary.tree.pdf", useDingbats = F, height = 5)
		for(i in 1:length(empstats)){
		      sdstat <- sd(statsmat[2:nrow(statsmat), i])
		      hist(statsmat[2:nrow(statsmat), i], xlim = c(min(statsmat[, i], na.rm = T) - sdstat, max(statsmat[, i], na.rm = T) + sdstat), xlab = statlabels[i], ylab = "Frequency of predictive simulations", main = "")
		      abline(v = empstats[i], col = "red", lwd = 3)
		}
		dev.off()
		
		
	}
	
	if("allqp" %in% whatToOutput){
	write.csv(geneResults[[1]], file = "full.results.csv")

if("pvals" %in% whatToOutput){
	rownames(geneResults[[1]]) <- geneResults[[1]][,1]
	restab <- geneResults[[1]][,-grep("ID|sCF|sDF|sN|gCF|gDF|gN", colnames(geneResults[[1]]))]
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