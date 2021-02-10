processTime <- proc.time()[3]
source("otherScripts/R/clean.gene.R")
source("otherScripts/R/test.saturation.R")
source("otherScripts/R/runIQtree.R")
source("testStatistics/get.entropy.test.R")
load("testStatistics/thresholdFunctions.Rdata")

initial.dir <- getwd()

print("Functions required were loaded successfully")

outs <- list()

locilengths <- vector()

if(input$Ncores > 1) parallelise <- T else parallelise <- F

setwd(input$outputFolder)

print("Output folder was identified successfully")

if(!input$overwrite && (file.exists("saturation.test.phylomad") | file.exists("saturation.plots.pdf") | file.exists("multigene.saturation.plot.pdf"))){
       stop("Exisitng files will not be overwritten. Aborting.")
}

whatToOutput <- unlist(input$whatToOutput)

if(input$Ncores > 1) parallelise <- T else parallelise <- F

if("satPlots" %in% whatToOutput | "multiSatPlots" %in% whatToOutput) plotdat <- T else plotdat <- F

geneResults <- try(test.saturation(loci = as.character(input$dataPath[, 4]), iqtreePath = iqtreePath, para = parallelise, ncore = input$Ncores, clean = input$dataTreatment, stats = input$saturationStats, plotdat = plotdat, linmods = funclist))
failedLoci <- c(which(sapply(geneResults, class) == "try-error"), which(sapply(geneResults, is.null)))

locinames <- as.character(input$dataPath[, 1])

if(length(failedLoci) > 0){
	geneResults <- geneResults[-failedLoci]
	locinames <- locinames[-failedLoci]
	print(paste("Analysis of locus", failedLoci, "failed"))
}

if(length(geneResults) > 1){
	restabs <- lapply(geneResults, function(x) x[[1]])
	restab <- do.call(rbind, restabs)
} else {
        restab <- geneResults[[1]][[1]]
}

if(input$dataTreatment == "codonpos") rownames(restab) <- as.character(sapply(locinames, function(x) paste(x, c("pos1and2", "pos3"), sep = "_"))) else rownames(restab) <- locinames
colnames(restab) <- gsub("enth", "Entropy", colnames(restab))
colnames(restab) <- gsub("enthvar", "EntropyOnVarSites", colnames(restab))
for(i in grep("_TPR_", colnames(restab))) restab[which(restab[, i] > 1), i] <- 1
for(i in grep("_FPR_", colnames(restab))) restab[which(restab[, i] < 0), i] <- 0

if("tsat" %in% whatToOutput) write.csv(restab, file = "saturation.test.results.csv")

if("satPlots" %in% whatToOutput){
	pdf("saturation.plots.pdf", height = 5, width = 5, useDingbats = F)
	for(i in 1:length(locinames)){
		if(input$dataTreatment == "codonpos"){
			plot(geneResults[[i]][[2]][[2]][[1]] ~ geneResults[[i]][[2]][[2]][[2]], pch = 20, col = "purple", ylab = "Uncorrected pairwise genetic distances", xlab = "Pairwise distances including the ratio\nof transitions to transversions (Tamura and Nei 1993)", main = paste0("Saturation plot for\n", locinames[i]))
			lm12 <- lm(geneResults[[i]][[2]][[2]][[1]] ~ 0 + geneResults[[i]][[2]][[2]][[2]])
			abline(lm12, lwd=2, col = "purple")
			points(geneResults[[i]][[2]][[1]][[1]] ~ geneResults[[i]][[2]][[1]][[2]], pch = 20, col = "red")
			lm3 <- lm(geneResults[[i]][[2]][[1]][[1]] ~ 0 + geneResults[[i]][[2]][[1]][[2]])
			abline(lm3, lwd=2, col = "red")
			cor12 <- round(cor.test(geneResults[[i]][[2]][[1]][[1]], geneResults[[i]][[2]][[1]][[2]])$estimate, 3)
			cor3 <- round(cor.test(geneResults[[i]][[2]][[2]][[1]], geneResults[[i]][[2]][[2]][[2]])$estimate, 3)
			legend("bottomright", legend = c(paste0("Pos 1+2, cor = ", cor12, ", slope = ", round(coef(lm12)[1], 3)), paste0("Pos 3, cor = ", cor3, ", slope = ", round(coef(lm3)[1], 3))), lty = 1, lwd = 2, col = c("red", "purple"), cex = 0.7)
		} else {
			plot(geneResults[[i]][[2]][[1]][[1]] ~ geneResults[[i]][[2]][[1]][[2]], pch = 20, col = "red", ylab = "Uncorrected pairwise genetic distances", xlab = "Pairwise distances including the ratio\nof transitions to transversions (Tamura and Nei 1993)", main = paste0("Saturation plot for\n", locinames[i]))
			lmdat <- lm(geneResults[[i]][[2]][[1]][[1]] ~ 0 + geneResults[[i]][[2]][[1]][[2]])
			abline(lmdat, lwd=2, col = "red")
			corlocus <- round(cor.test(geneResults[[i]][[2]][[1]][[1]], geneResults[[i]][[2]][[1]][[2]])$estimate, 3)
			legend("bottomright", legend = paste0("cor = ", corlocus, ", slope = ", round(coef(lmdat)[1], 3)), cex = 0.7)
		}
		abline(0, 1, lty=2)
	}
	dev.off()
}

if("multiSatPlots" %in% whatToOutput){
	pdf("multilocus.plots.pdf", height = 5, width = 5, useDingbats = F)
	risklegend <- c("Low", "Medium", "High")
	if(input$dataTreatment == "codonpos"){
	      colsplot <- c("red", "purple")
	      risklegcol <- c(rep("black", 5), "red", "purple")
	} else {
	      colsplot <- rep("black", nrow(restab))
	      risklegcol <- c("black", "blue", "black", "red")
	}
	if("t_Entropy" %in% colnames(restab)){
	      resprops <- c(length(which(restab[, "Risk_Entropy"] == "LOW")), length(which(restab[, "Risk_Entropy"] == "MEDIUM")), length(which(restab[, "Risk_Entropy"] == "HIGH")))/nrow(restab)
	      risklegendent <- paste0(risklegend, " (", as.character(round(resprops, 2)), ")")
	      if(input$dataTreatment == "codonpos"){
	      		risklegendent <- c("Risk (prop. loci)", risklegendent, "Codon position", "1+2", "3")
	      } else {
	      		risklegendent <- c("Risk (prop. loci)", risklegendent)
			colsplot[which(restab[, "Risk_Entropy"] == "HIGH")] <- "red"
			colsplot[which(restab[, "Risk_Entropy"] == "LOW")] <- "blue"
	      }
	      pchsplot <- restab[,"Risk_Entropy"]
	      pchsplot[which(pchsplot == "LOW")] <- 16
	      pchsplot[which(pchsplot == "MEDIUM")] <- 17
	      pchsplot[which(pchsplot == "HIGH")] <- 18
	      plot(restab[, "N_sites"], restab[, "t_Entropy"], col = colsplot, pch = as.numeric(pchsplot), main = "Multilocus saturation test\nEntropy statistic", xlab = "N sites", ylab = "t (entropy statistic)")
	      legend("topleft", pch = c(NA, 16:18, NA, 16, 16), legend = risklegendent, col = risklegcol)
	}
	if("t_CI" %in% colnames(restab)){
	      resprops <- c(length(which(restab[, "Risk_EntropyVarSites"] == "LOW")), length(which(restab[, "Risk_EntropyVarSites"] == "MEDIUM")), length(which(restab[, "Risk_EntropyVarSites"] == "HIGH")))/nrow(restab)
	      risklegendci <- paste0(risklegend, " (", as.character(round(resprops, 2)), ")")
              if(input$dataTreatment == "codonpos"){
                        risklegendci <- c("Risk (prop. loci)", risklegendci, "Codon position", "1+2", "3")
              } else {
			risklegendci <- c("Risk (prop. loci)", risklegendci)
	      		colsplot <- rep("black", nrow(restab))
			colsplot[which(restab[, "Risk_EntropyVarSites"] == "HIGH")] <- "red"
			colsplot[which(restab[, "Risk_EntropyVarSites"] == "LOW")] <- "blue"
              }
              pchsplot <- restab[, "Risk_EntropyVarSites"]
              pchsplot[which(pchsplot == "LOW")] <- 16
              pchsplot[which(pchsplot == "MEDIUM")] <- 17
              pchsplot[which(pchsplot == "HIGH")] <- 18
              plot(restab[, "N_sites"], restab[, "t_EntropyVarSites"], col = colsplot, pch = as.numeric(pchsplot), main = "Multilocus saturation test\nEntropy on var sites", xlab = "N sites", ylab = "t (entropy on var sites)")
              legend("topleft", pch = c(NA, 16:18, NA, 16, 16), legend = risklegendentvar, col = risklegcol)
        }
	dev.off()
}

setwd(initial.dir)

print(paste0("Assessment completed successfully in ", round(proc.time()[3] - processTime, 3), " seconds."))