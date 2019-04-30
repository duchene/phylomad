processTime <- proc.time()[3]

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
       system(paste0("mkdir ", as.character(input$dataPath[j, 1]), ".phylomad.phylosig"))
       setwd(paste0(as.character(input$dataPath[j, 1]), ".phylomad.phylosig"))
}

whatToOutput <- unlist(input$whatToOutput)

geneResults <- try(test.phylosignal(sdata = analysisdata, format = if(input$testType == "locus") "bin" else dataFormat, testType = input$testType, aadata = input$dataType, model = model, iqtreePath = iqtreePath, astralPath = astralPath, Nsims = input$Nsims, testStats = selectedStats, returnEstPhylo = "phyloempres" %in% whatToOutput, returnSimulations = "simdat" %in% whatToOutput))
if(class(geneResults) == "try-error"){
       setwd("..")
       system(paste0("rm -r ", as.character(input$dataPath[j, 1]), ".phylomad.phylosig"))
       print(paste0("Assessment of ", as.character(input$dataPath[j, 1]), " failed"))
       next
}

rownames(geneResults[[1]]) <- geneResults[[1]][,1]
colnames(geneResults[[1]]) <- gsub("Length", "br.length", colnames(geneResults[[1]]))
colnames(geneResults[[1]]) <- gsub("Label", "br.support", colnames(geneResults[[1]]))
geneResults[[1]] <- round(as.data.frame(apply(geneResults[[1]], 2, as.numeric)), 3)
geneResults[[1]] <- rbind(geneResults[[1]], round(colMeans(geneResults[[1]], na.rm = T), 3))
rownames(geneResults[[1]])[nrow(geneResults[[1]])] <- "mean"

if("allqp" %in% whatToOutput){
	write.csv(t(geneResults[[1]]), file = "full.results.csv")
}

if("pvals" %in% whatToOutput){
	restab <- geneResults[[1]][,-grep("ID|sDF|sN|gDF|gN", colnames(geneResults[[1]]))]
	write.csv(t(restab), file = "results.summary.csv")
}

if("testPlots" %in% whatToOutput){
	histplotdat <- as.numeric(geneResults[[1]]["mean",])
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
	        geneResults[[1]] <- as.matrix(geneResults[[1]])
		pdf("tests.density.plots.pdf", useDingbats = F, height = 4, width = 15)
		par(mfrow = c(1, 3))
		for(i in 1:length(selectedStats)){
		      statdat <- histplotdat[grep(selectedStats[i], names(histplotdat))]
		      if(is.na(statdat[1]) || is.nan(statdat[1]) || any(is.infinite(statdat)) || all(is.na(statdat)) || all(is.nan(statdat))){
		      	print(paste0("Plots of ", statlabels[i], "cannot be created."))
		      	next
		      }
		      
		      # Plot of all simulated branch values vs. all empirical branch values
		      statsim <- geneResults[[1]][-nrow(geneResults[[1]]), grep(selectedStats[i], colnames(geneResults[[1]]))]
		      statsim <- as.numeric(statsim[,2:(input$Nsims+1)])
		      statemp <- as.numeric(geneResults[[1]][-nrow(geneResults[[1]]), if(selectedStats[i] == "CF") 2 else paste0("emp.", selectedStats[i])])
		      statsimdens <- try(density(statsim, from = min(statsim, na.rm = T), to = max(statsim, na.rm = T), na.rm = T))
		      statempdens <- try(density(statemp, from = min(statemp, na.rm = T), to = max(statemp, na.rm = T), na.rm = T))
		      if(class(statsimdens) != "try-error" & class(statempdens) != "try-error"){
		      		plot(statsimdens, xlim = range(c(statsim, statemp), na.rm = T), ylim = c(0, max(c(statsimdens$y, statempdens$y), na.rm = T)), xlab = statlabels[i], main = "")
		      		lines(statempdens, col = "red")
		      		legend("topright", legend = c("Empirical branches", "Simlulated branches"), col = c("red", "black"), lty = 1)
		      } else { 
				frame() 
		      }
		      
		      # Plot of branch-wise p-values
		      statpvals <- as.numeric(geneResults[[1]][-nrow(geneResults[[1]]), paste0(selectedStats[i], ".p.value")])
		      pdens <- try(density(statpvals, from = min(statpvals, na.rm = T), to = max(statpvals, na.rm = T), na.rm = T))
		      if(class(pdens) != "try-error") plot(pdens, main = "", xlab = paste0("Branch P-values\n(", statlabels[i], ")")) else frame()
		      
		      # Plot of test of mean across branches
		      statmeandens <- try(density(statdat[-1], from = min(statdat[-1], na.rm = T), to = max(statdat[-1], na.rm = T), na.rm = T))
		      if(class(statmeandens) != "try-error"){
		      		plot(statmeandens, main = "", xlim = range(statdat, na.rm = T), xlab = paste0("Mean ", statlabels[i], " across branches"))
				abline(v = statdat[1], lwd = 2, col = "red")
				legend("topright", legend = c("Empirical", "Simulated"), lty = 1, col = c("red", "black"))
		      } else {
				frame()
		      }
		      
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

		      plotBranchbyTrait(geneResults[[2]], brpvalue, mode = "edges", palette = "rainbow", type = "unrooted", legend = F)
		      plotBranchbyTrait(tr, brpvalue, mode = "edges", palette = "rainbow", type = "unrooted", title = paste0(statlabels[i], "\nP-value\n"))
		      plotBranchbyTrait(geneResults[[2]], brsdpd, mode = "edges", palette = "rainbow", type = "unrooted", legend = F)
                      plotBranchbyTrait(tr, brsdpd, mode = "edges", palette = "rainbow", type = "unrooted", title = paste0(statlabels[i], "\nz-score\n"))
		}
		dev.off()
		
		#finish up and add custom legend? Add two ternary plots one with all the simulated branches, and one of the mean cfs for each data set, put both in one page with grid.arrange
		pdf("ternary.plots.pdf", height = 4, width = 8)
		#geneResults[[1]] <- as.data.frame(geneResults[[1]])
		if(input$testType == "locus") ttyp <- "s" else ttyp <- "g"
		terndat <- geneResults[[1]][-nrow(geneResults[[1]]), paste0(ttyp, c("CF", "DF1", "DF2"))]
		ternmeandat <- colMeans(terndat)
		nbranch <- nrow(terndat)
		nsimbranch <- nbranch * input$Nsims
		for(i in 1:input$Nsims){
		      terndat <- rbind(geneResults[[1]][-nrow(geneResults[[1]]),c(paste0(ttyp, "CF.sim.", i), paste0(ttyp, "DF1.sim.", i), paste0(ttyp, "DF2.sim.", i))], terndat)
		      ternmeandat <- rbind(geneResults[[1]][nrow(geneResults[[1]]),c(paste0(ttyp, "CF.sim.", i), paste0(ttyp, "DF1.sim.", i), paste0(ttyp, "DF2.sim.", i))], ternmeandat)
		}
		colnames(terndat) <- colnames(ternmeandat) <- c("CF", "DF1", "DF2")
		
		terndat <- as.data.frame(terndat)
		#terndat$cols <- c(rep(2, nsimbranch), rep(1, nbranch))
		terndat$cols <- as.factor(rep(1:nbranch, input$Nsims+1))
		ternallbrs <- ggtern(data=terndat, aes(x = DF1, y = CF, z = DF2),aes(x, y, z)) + geom_point(aes(fill = cols), alpha = c(rep(0.25, nsimbranch), rep(1, nbranch)), stroke = 0, size = 2, shape = c(rep(21, nsimbranch), rep(23, nbranch))) + theme(legend.position = "none") + ggtitle("All simulated and empirical\nbranches")
		
		ternmeandat <- as.data.frame(ternmeandat)
		ternmeandat$cols <- c(rep(2,input$Nsims), 1)
		ternmean <- ggtern(data=ternmeandat,aes(x=DF1,y=CF,z=DF2),aes(x,y,z)) + geom_point(aes(fill= cols), alpha = c(rep(0.5,input$Nsims), 1), stroke=0,size=2,shape=c(rep(21,input$Nsims), 23)) + theme(legend.position = "none") + ggtitle("Mean across each simulation\nand empirical tree")
		
		grid.arrange(ternallbrs, ternmean, nrow = 1)
		dev.off()
		
	}
}

setwd("..")

}

setwd(initial.dir)

print(paste0("Assessment completed successfully in ", round(proc.time()[3] - processTime, 3), " seconds."))