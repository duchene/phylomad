processTime <- proc.time()[3]

initial.dir <- getwd()

source("otherScripts/R/test.phylosignal.R")
source("otherScripts/R/runIQtree.R")
source("otherScripts/R/clean.gene.R")
source("otherScripts/R/get.model.R")
source("otherScripts/R/sample.quartet.R")
source("testStatistics/get.phylosignal.metrics.R")
setwd("otherScripts/phybase/")
for(i in dir()) source(i)
setwd(initial.dir)
print("Functions required were loaded successfully")

selectedStats <- unlist(input$testStats)
outs <- list()
locilengths <- vector()

if(input$testType == "tree"){
	firstLine <- readLines(as.character(input$sptreePath[1, 4]), n = 1)
	if(grepl("[#]NEXUS|[#]nexus", firstLine)) treesFormat <- "nexus" else if(grepl("[(]", firstLine)) treesFormat <- "newick"
	if(treesFormat == "newick"){ sptreraw <- read.tree(as.character(input$sptreePath[1, 4])) } else if(treesFormat == "nexus"){ sptreraw <- read.nexus(as.character(input$sptreePath[1, 4])) }
	if(is.rooted(sptreraw)) sptre <- unroot(sptreraw) else sptre <- sptreraw
	# List of edges, each with with a list of Nsims quartets
	allsampquarts <- lapply(which(!sptre$edge[,2] %in% 1:Ntip(sptre)), function(x) lapply(1:input$Nsims, function(y) sample.quartet(sptre, x)))
	# Matrices where columns are edges and rows are replicates
        medmat <- matrix(NA, ncol = length(allsampquarts), nrow = input$Nsims * nrow(input$dataPath))
        pvalmat <- matrix(NA, ncol = length(allsampquarts), nrow = input$Nsims * nrow(input$dataPath))
        presabs <- matrix(NA, ncol = length(allsampquarts), nrow = input$Nsims * nrow(input$dataPath))
        colnames(medmat) <- colnames(pvalmat) <- colnames(presabs) <- names(allsampquarts) <- which(!sptre$edge[,2] %in% 1:Ntip(sptre))
        rownames(medmat) <- rownames(pvalmat) <- rownames(presabs) <- as.character(sapply(1:nrow(input$dataPath), function(x) paste0("rep.", 1:input$Nsims, ".locus.", x)))
}

for(j in 1:nrow(input$dataPath)){
      firstLine <- readLines(as.character(input$dataPath[j, 4]), n = 1)
      if(grepl("[>]", firstLine)) dataFormat <- "fasta" else if(grepl("[#]NEXUS|[#]nexus", firstLine)) dataFormat <- "nexus" else if(grepl("[(]", firstLine)) dataFormat <- "newick" else dataFormat <- "phylip"

      #if(input$testType %in% c("locus", "tree") & !input$dataType) model <- get.model(as.character(input$dataPath[j, 4]), format = dataFormat) else if(input$testType %in% c("locus", "tree") & input$dataType) model <- "WAG+G"
      
      model <- "GTR+G"

      print("Model to be assessed was identified")

      if(input$testType %in% c("locus", "tree")) analysisdata <- clean.gene(sdata = as.character(input$dataPath[j, 4]), format = dataFormat, aadata = input$dataType, clean = F) else analysisdata <- as.character(input$dataPath[j, 4])

      print("Locus was cleaned successfully")

      setwd(input$outputFolder)

      print("Output folder was identified successfully")

      whatToOutput <- unlist(input$whatToOutput)

      if(input$testType == "locus"){
      
      if(!input$overwrite && file.exists(paste0(as.character(input$dataPath[j, 1]), ".phylomad.phylosig"))){
      	stop("Exisitng files will not be overwritten. Aborting.")
      } else {
       	system(paste0("mkdir ", as.character(input$dataPath[j, 1]), ".phylomad.phylosig"))
       	setwd(paste0(as.character(input$dataPath[j, 1]), ".phylomad.phylosig"))
      }

      geneResults <- try(test.phylosignal(sdata = analysisdata, format = if(input$testType == "locus") "bin" else dataFormat, testType = input$testType, aadata = input$dataType, model = model, iqtreePath = iqtreePath, astralPath = astralPath, Nsims = input$Nsims, testStats = selectedStats, returnSimulations = "simdat" %in% whatToOutput))
      if(class(geneResults) == "try-error"){
       	setwd("..")
       	system(paste0("rm -r ", as.character(input$dataPath[j, 1]), ".phylomad.phylosig"))
       	print(paste0("Assessment of ", as.character(input$dataPath[j, 1]), " failed"))
       	next
      }

      rownames(geneResults[[1]]) <- geneResults[[1]][,1]
      colnames(geneResults[[1]]) <- gsub("Length", "br.length", colnames(geneResults[[1]]))
      colnames(geneResults[[1]]) <- gsub("Label", "br.support", colnames(geneResults[[1]]))
      for(y in 2:ncol(geneResults[[1]])) geneResults[[1]][,y] <- round(as.numeric(geneResults[[1]][,y]), 3)
      geneResults[[1]] <- rbind(geneResults[[1]], round(colMeans(geneResults[[1]], na.rm = T), 3))
      rownames(geneResults[[1]])[nrow(geneResults[[1]])] <- "mean"

      ## Annotate empirical tree with bootstrap support ranges and p-valuea
      if(input$testType == "locus" & "phyloempres" %in% whatToOutput){
	nexhead	<- c("#NEXUS", "begin trees;")
	antre <- geneResults[[2]]
	antre$node.label <- vector()
	for(y in 1:(nrow(geneResults[[1]] - 1))){
	      annot <- ""
	      for(z in 1:length(selectedStats)){
	      	    fullstat <- if(selectedStats[z] == "CF") as.numeric(geneResults[[1]][y, paste0("sCF.sim.", 1:input$Nsims)]) else as.numeric(geneResults[[1]][y, paste0("sim.", 1:input$Nsims, ".", selectedStats[z])])
	      	    statmed <- round(median(fullstat, na.rm = T), 2)
		    statran <- round(quantile(fullstat, probs = c(0.025, 0.975), na.rm = T), 2)
	      	    annot <- c(annot, paste0(selectedStats[z], '=\"', statmed, ' [', statran[1], ", ", statran[2], '] ovcon.p=', geneResults[[1]][y, paste0(selectedStats[z], ".p.value")], '\"'))
	      }
	      annot <- paste0("<&!", paste0(annot, collapse = ","), ">")
	      antre$node.label[antre$edge[which(geneResults[[3]]$edge.length == geneResults[[1]][y, 1]), 2]-Ntip(antre)] <- annot
	      
	}
	antre <- gsub("[_]", " ", write.tree(antre))
	antre <- gsub("[-]", ",", antre)
	antre <- gsub("[[]", "(", antre)
	antre <- gsub("[]]", ")", antre)
	antre <- gsub("[<]", "[", antre)
        antre <- gsub("[>]", "]", antre)
	antre <- c(nexhead, paste0("tree tree_1 = [&R] ", antre), "end;")
	writeLines(antre, con = "estimated.tree.supports.tre")
}

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
		      	print(paste0("Plots of ", statlabels[i], " cannot be created."))
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
		
		pdf("tests.summary.tree.pdf", useDingbats = F, height = if(length(geneResults[[2]]$edge.length) < 50) 5 else 10, width = if(length(geneResults[[2]]$edge.length) < 50) 10 else 20)
		for(i in 1:length(selectedStats)){
		      par(mfrow = c(1,2), mar = c(5.1, 4.1, 4.1, 2.1))
		      tr <- geneResults[[2]]
		      tr$edge.length <- rep(1, length(geneResults[[2]]$edge.length))
		      brpvalue <- geneResults[[1]][, paste0(selectedStats[i], ".p.value")]
		      brpvalue[which(brpvalue > 0.01)] <- 1
		      brpvalue[which(brpvalue <= 0.01)] <- 2
		      brsdpd <- round(geneResults[[1]][, paste0(selectedStats[i], ".sdpd")], 1)
		      names(brpvalue) <- names(brsdpd) <- geneResults[[1]][,1]
		      
		      brpvalue <- brpvalue[as.character(geneResults[[3]]$edge.length)]
		      brsdpd <- brsdpd[as.character(geneResults[[3]]$edge.length)]
		      brpvalue[is.na(brpvalue)] <- 1 # mean(brpvalue, na.rm = T)
		      brsdpd[is.na(brsdpd)] <- mean(brsdpd, na.rm = T)
		      if(all(is.na(brpvalue)) || all(is.nan(brpvalue)) || any(is.infinite(brpvalue)) || all(is.na(brsdpd)) || all(is.nan(brsdpd)) || any(is.infinite(brsdpd))){
		      	print(paste0("Tree depicting ", statlabels[i], " statistic cannot be ploted."))
			next
		      }

		      #plotBranchbyTrait(geneResults[[2]], brpvalue, mode = "edges", palette = "rainbow", type = "unrooted", legend = F)
		      #plotBranchbyTrait(tr, brpvalue, mode = "edges", palette = "rainbow", type = "unrooted", title = paste0(statlabels[i], "\nP-value\n"))
		      plot(geneResults[[2]], edge.color = brpvalue, type = "unrooted", main = paste0(statlabels[i], ". In red P-values < 0.01\n(includes branch lengths)"))
		      plot(tr, edge.color = brpvalue, type = "unrooted", main = "\n(excludes branch lengths)")
		      brpvalcols <- brpvalue
		      brpvalcols[brpvalcols == 1] <- "white"
		      brpvalcols[brpvalcols == "2"] <- "red"
		      edgelabels(names(brpvalue), frame = "circle", bg = brpvalcols, cex = 0.7)
		      
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

} else if(input$testType == "tree"){

       if(!input$overwrite && file.exists(paste0(as.character(input$sptreePath[1, 1]), ".phylomad.supports.tre"))) stop("Exisitng files will not be overwritten. Aborting.")

	if(input$Ncores == 1){

	for(x in 1:length(allsampquarts)){
	      allqtax <- sapply(allsampquarts[[x]], function(y) all(y$tip.label %in% rownames(analysisdata)))
              if(!all(allqtax)){
		print("Quartet missing from tree")
                next
              }
	      for(z in which(allqtax)){
	      	    quartloc <- analysisdata[allsampquarts[[x]][[z]]$tip.label,]
	      	    geneResult <- try(test.phylosignal(sdata = quartloc, format = "bin", testType = "locus", aadata = input$dataType, model = model, iqtreePath = iqtreePath, astralPath = astralPath, Nsims = input$Nsims, testStats = "CF", returnSimulations = F))
		    fullstat <- as.numeric(geneResult[[1]][1, paste0("sCF.sim.", 1:input$Nsims)])
            medmat[z + ((j-1) * input$Nsims), x] <- round(median(fullstat, na.rm = T), 3)
		    pvalmat[z + ((j-1) * input$Nsims), x] <- round(geneResult[[1]][1, "CF.p.value"], 3)
		    presabs[z + ((j-1) * input$Nsims), x] <- (RF.dist(allsampquarts[[x]][[z]], geneResult[[2]]) != 0)
	      }
	      pvalmat[which(allqtax) + ((j-1) * input$Nsims), x] <- round(p.adjust(pvalmat[which(allqtax) + ((j-1) * input$Nsims), x], method = "fdr"), 3)   
	}
	
	} else {
	
		### START PARALLEL COMPUTING
	   print('Parallel computing started')
	   system(paste0("mkdir processing.locus.", j))
	   setwd(paste0("processing.locus.", j))
	   require(foreach)
	   require(doParallel)
	   dataType <- input$dataType
	   Nsims <- input$Nsims
	   runSim <- function(x){
	   		system(paste0("mkdir branchfolder.", x))
	   		setwd(paste0("branchfolder.", x))
	   		allqtax <- sapply(allsampquarts[[x]], function(y) all(y$tip.label %in% rownames(analysisdata)))
              if(!all(allqtax)){
		print("Quartet missing from tree")
                next
              }
	      medmatdat <- matrix(NA, ncol = 1, nrow = length(allsampquarts[[x]]))
	      pvalmatdat <- matrix(NA, ncol = 1, nrow = length(allsampquarts[[x]]))
	      presabsdat <- matrix(NA, ncol = 1, nrow = length(allsampquarts[[x]]))	      
	      for(z in which(allqtax)){
	      	    quartloc <- analysisdata[allsampquarts[[x]][[z]]$tip.label,]
	      	    geneResult <- try(test.phylosignal(sdata = quartloc, format = "bin", testType = "locus", aadata = dataType, model = model, iqtreePath = iqtreePath, astralPath = astralPath, Nsims = Nsims, testStats = "CF", returnSimulations = F))
		    fullstat <- as.numeric(geneResult[[1]][1, paste0("sCF.sim.", 1:Nsims)])
            medmatdat[z, 1] <- round(median(fullstat, na.rm = T), 3)
		    pvalmatdat[z, 1] <- round(geneResult[[1]][1, "CF.p.value"], 3)
		    presabsdat[z, 1] <- (RF.dist(allsampquarts[[x]][[z]], geneResult[[2]]) != 0)
	      }
	      pvalmatdat[which(allqtax), 1] <- round(p.adjust(pvalmatdat[which(allqtax), 1], method = "fdr"), 3)
	      setwd("..")
	      system(paste0("rm -r branchfolder.", x))
	      res <- list(medmatdat, pvalmatdat, presabsdat)
	      return(res)
	   }
	   
	   cl <- makeCluster(input$Ncores)
	   registerDoParallel(cl)
	   simReps <- foreach(x = 1:length(allsampquarts), .packages = c('phangorn', 'ape'), .export = c('test.phylosignal', 'runIQtree')) %dopar% runSim(x)
	   for(x in 1:length(allsampquarts)){
	   		medmat[(1:input$Nsims) + ((j-1) * input$Nsims), x] <- simReps[[x]][[1]]
		    pvalmat[(1:input$Nsims) + ((j-1) * input$Nsims), x] <- simReps[[x]][[2]]
		    presabs[(1:input$Nsims) + ((j-1) * input$Nsims), x] <- simReps[[x]][[3]]
	   } 
	   stopCluster(cl)
	   setwd("..")
	   system(paste0("rm -r processing.locus.", j))
	   print("Parallel computing ended successfully")
	   ### END PARALLEL COMPUTING
	   
	}
	
	
	if("pvals" %in% whatToOutput){
	    write.csv(medmat, file = paste0(as.character(input$sptreePath[1, 1]), ".phylomad.median.CFs.csv"))
	    write.csv(pvalmat, file = paste0(as.character(input$sptreePath[1, 1]), ".phylomad.median.pvals.csv"))
	    write.csv(presabs, file = paste0(as.character(input$sptreePath[1, 1]), ".phylomad.median.presence.csv"))
		}
	if("phyloempres" %in% whatToOutput){
	nexhead <- c("#NEXUS", "begin trees;")
	annot <- sapply(1:length(allsampquarts), function(x){
	      medmed <- round(median(medmat[,x], na.rm = T), 2)
	      ranmed <- round(quantile(medmat[,x], probs = c(0.025, 0.975), na.rm = T), 2)
	      medpval <- round(median(pvalmat[,x], na.rm = T), 2)
	      ranpval <- round(quantile(pvalmat[,x], probs = c(0.025, 0.975), na.rm = T), 2)
	      medabs <- round(mean(presabs[,x], na.rm = T), 2)
	      ranabs <- round(quantile(presabs[,x], probs = c(0.025, 0.975), na.rm = T), 2)
	      ann <- paste0('<&!,CF=\"', medmed, ' [', ranmed[1], ", ", ranmed[2], ']\",OverconP=\"', medpval, ' [', ranpval[1], ", ", ranpval[2], ']\",Absence=\"', medabs, ' [', ranabs[1], ", ", ranabs[2], ']\">')
	      return(ann)
	})
	antre <- sptre
	antre$node.label <- c(NA, annot)
	antre <- gsub("[_]", " ", write.tree(antre))
        antre <- gsub("[-]", ",", antre)
        antre <- gsub("[[]", "(", antre)
        antre <- gsub("[]]", ")", antre)
        antre <- gsub("[<]", "[", antre)
        antre <- gsub("[>]", "]", antre)
        antre <- c(nexhead, paste0("tree tree_1 = [&R] ", antre), "end;")
        writeLines(antre, con = paste0(as.character(input$sptreePath[1, 1]), ".phylomad.support.tre"))
	}
	
}
}

setwd(initial.dir)

print(paste0("Assessment completed successfully in ", round(proc.time()[3] - processTime, 3), " seconds."))