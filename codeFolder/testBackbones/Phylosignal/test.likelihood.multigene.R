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

firstLine <- readLines(as.character(input$dataPath[j, 4]), n = 1)
if(grepl("[>]", firstLine)) dataFormat <- "fasta" else if(grepl("[#]NEXUS|[#]nexus", firstLine)) dataFormat <- "nexus" else if(grepl("[(]", firstLine)) dataFormat <- "newick" else dataFormat <- "phylip"

if(input$testType == "locus" & !input$dataType) model <- get.model(as.character(input$dataPath[j, 4]), format = dataFormat) else if(input$testType == "locus" & input$dataType) model <- "WAG+G"

print("Model to be assessed was identified")

if(input$testType == "locus") analysisdata <- clean.gene(sdata = as.character(input$dataPath[j, 4]), format = dataFormat, aadata = input$dataType, clean = F) else analysisdata <- as.character(input$dataPath[j, 4])

print("Locus was cleaned successfully")

setwd(input$outputFolder)

print("Output folder was identified successfully")

if(!input$overwrite && file.exists(paste0(as.character(input$dataPath[j, 1]), ".phylomad"))){
       stop("Exisitng files will not be overwritten. Aborting.")
} else {
       system(paste0("mkdir ", as.character(input$dataPath[j, 1]), ".phylomad"))
       setwd(paste0(as.character(input$dataPath[j, 1]), ".phylomad"))
}

whatToOutput <- unlist(input$whatToOutput)

geneResults <- test.phylosignal(sdata = analysisdata, format = if(input$testType == "locus") "bin" else dataFormat, testType = input$testType, aadata = input$dataType, model = model, iqtreePath = iqtreePath, astralPath = astralPath, Nsims = input$Nsims, testStats = selectedStats, returnEstPhylo = "phyloempres" %in% whatToOutput, returnSimulations = "simdat" %in% whatToOutput)

rownames(geneResults[[1]]) <- geneResults[[1]][,1]
colnames(geneResults[[1]]) <- gsub("Label", "br.length", colnames(geneResults[[1]]))
geneResults[[1]] <- round(as.data.frame(apply(geneResults[[1]], 2, as.numeric)), 3)

	if("allqp" %in% whatToOutput) write.csv(t(geneResults[[1]]), file = "full.results.csv")

if("pvals" %in% whatToOutput){
	restab <- geneResults[[1]][,-grep("ID|sDF|sN|gDF|gN", colnames(geneResults[[1]]))]
	restab <- rbind(restab, round(colMeans(restab, na.rm = T), 3))
	rownames(restab)[nrow(restab)] <- "mean"
	write.csv(t(restab), file = "results.summary.csv")
}

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
		      if(is.na(statdat[1]) || is.nan(statdat[1]) || any(is.infinite(statdat)) || all(is.na(statdat)) || all(is.nan(statdat))){
		      	print(paste0("Histogram of ", statlabels[i], "cannot be plotted."))
		      	next
		      }
		      sdstat <- sd(statdat[2:length(statdat)], na.rm = T)
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
		      if(all(is.na(brpvalue)) || all(is.nan(brpvalue)) || any(is.infinite(brpvalue)) || all(is.na(brsdpd)) || all(is.nan(brsdpd)) || any(is.infinite(brsdpd))){
		      	print(paste0("Tree depicting ", statlabels[i], " statistic cannot be ploted."))
			next
		      }

		      plotBranchbyTrait(geneResults[[2]], brpvalue, mode = "edges", palette = "heat.colors", type = "fan", legend = F)
		      plotBranchbyTrait(tr, brpvalue, mode = "edges", palette = "heat.colors", type = "fan", title = paste0(statlabels[i], "\nP-value\n"))
		      plotBranchbyTrait(geneResults[[2]], brsdpd, mode = "edges", palette = "heat.colors", type = "fan", legend = F)
                      plotBranchbyTrait(tr, brsdpd, mode = "edges", palette = "heat.colors", type = "fan", title = paste0(statlabels[i], "\nz-score\n"))
		}
		dev.off()
		
		
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