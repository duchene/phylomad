source("otherScripts/R/clean.gene.R")
source("otherScripts/R/test.saturation.R")
source("otherScripts/R/runPhyML.R")
source("testStatistics/get.entropy.test.R")
source("testStatistics/get.ci.test.R")
source("testStatistics/get.comp.test.R")
load("testStatistics/thresholdFunctions.Rdata")

initial.dir <- getwd()

print("Functions required were loaded successfully")

outs <- list()

locilengths <- vector()

print("Locus was cleaned successfully")

if(input$Ncores > 1) parallelise <- T else parallelise <- F

setwd(input$outputFolder)

print("Output folder was identified successfully")

if(!input$overwrite && (file.exists("saturation.test.phylomad") | file.exists("saturation.plots.pdf") | file.exists("multigene.saturation.plot.pdf"))){
       stop("Exisitng files will not be overwritten. Aborting.")
}

whatToOutput <- unlist(input$whatToOutput)

if(input$Ncores > 1) parallelise <- T else parallelise <- F

if("satPlots" %in% whatToOutput | "multiSatPlots" %in% whatToOutput) plotdat <- T else plotdat <- F

geneResults <- test.saturation(loci = as.character(input$dataPath[, 4]), format = input$dataFormat, phymlPath = phymlPath, para = parallelise, ncore = input$Ncores, clean = input$dataTreatment, stats = input$saturationStats, plotdat = plotdat, linmods = funclist)

#### Output missing

locinames <- as.character(input$dataPath[, 1])
	
if(length(geneResults) > 1){
	restabs <- lapply(geneResults, function(x) x[[1]])
	restab <- do.call(rbind, restabs)
} else {
        restab <- geneResults[[1]][[1]]
}
if(input$dataTreatment == "codonpos") rownames(restab) <- as.character(sapply(locinames, function(x) paste(x, c("pos1and2", "pos3"), sep = "_"))) else rownames(restab) <- locinames
colnames(restab) <- gsub("enth", "Entropy", colnames(restab))
colnames(restab) <- gsub("cith", "CI", colnames(restab))
colnames(restab) <- gsub("comth", "Compression", colnames(restab))

if("tsat" %in% whatToOutput){
	write.csv(restab, file = "saturation.test.results.csv")
}

if("satPlots" %in% whatToOutput){
	pdf("saturation.plots.pdf", height = 5, width = 5, useDingbats = F)
	for(i in 1:length(locinames)){
		if(input$dataTreatment == "codonpos"){
			plot(geneResults[[i]][[2]][[2]][[1]] ~ geneResults[[i]][[2]][[2]][[2]], pch = 20, col = "purple", ylab = "Uncorrected pairwise genetic distances", xlab = "Pairwise distances including the ratio\nof transitions to transversions (Tamura and Nei 1993)", main = paste0("Saturation plot for\n", locinames[i]))
			abline(lm(geneResults[[i]][[2]][[2]][[1]] ~ 0 + geneResults[[i]][[2]][[2]][[2]]), lwd=2, col = "purple")
			points(geneResults[[i]][[2]][[1]][[1]] ~ geneResults[[i]][[2]][[1]][[2]], pch = 20, col = "red")
			abline(lm(geneResults[[i]][[2]][[1]][[1]] ~ 0 + geneResults[[i]][[2]][[1]][[2]]), lwd=2, col = "red")
			cor12 <- round(cor.test(geneResults[[i]][[2]][[1]][[1]], geneResults[[i]][[2]][[1]][[2]])$estimate, 3)
			cor3 <- round(cor.test(geneResults[[i]][[2]][[2]][[1]], geneResults[[i]][[2]][[2]][[2]])$estimate, 3)
			legend("topleft", legend = c(paste0("1st and 2nd positions, cor = ", cor12), paste0("3rd position, cor = ", cor3)), lty = 1, lwd = 2, col = c("red", "purple"), cex = 0.7)
		} else {
			plot(geneResults[[i]][[2]][[1]][[1]] ~ geneResults[[i]][[2]][[1]][[2]], pch = 20, col = "red", ylab = "Uncorrected pairwise genetic distances", xlab = "Pairwise distances including the ratio\nof transitions to transversions (Tamura and Nei 1993)", main = paste0("Saturation plot for\n", locinames[i]))
			abline(lm(geneResults[[i]][[2]][[1]][[1]] ~ 0 + geneResults[[i]][[2]][[1]][[2]]), lwd=2, col = "red")
			corlocus <- round(cor.test(geneResults[[i]][[2]][[1]][[1]], geneResults[[i]][[2]][[1]][[2]])$estimate, 3)
			legend("topleft", legend = paste0("cor = ", corlocus), cex = 0.7)
		}
		abline(0,1, lty=2)
	}
	dev.off()
}

if("multiSatPlots" %in% whatToOutput){
	pdf("multilocus.plots.pdf", height = 7, width = 7, useDingbats = F)
	if(input$dataTreatment == "codonpos"){
	      colsplot <- c("red", "purple")
	      risklegend <- c("Risk", "Low", "Medium", "High", "Codon position", "1+2", "3")
	      risklegcol <- c(rep("black", 5), "red", "purple")
	} else {
	      colsplot <- rep("black", nrow(restab))
	      risklegend <- c("Risk", "Low", "Medium", "High")
	      risklegcol <- c("black", "blue", "black", "red")
	}
	if("t_Entropy" %in% colnames(restab)){
	      if(input$dataTreatment != "codonpos"){
	      		colsplot[which(restab[, "Risk_Entropy"] == "HIGH")] <- "red"
			colsplot[which(restab[, "Risk_Entropy"] == "LOW")] <- "blue"
	      }
	      pchsplot <- restab[,"Risk_Entropy"]
	      pchsplot[which(pchsplot == "LOW")] <- 16
	      pchsplot[which(pchsplot == "MEDIUM")] <- 17
	      pchsplot[which(pchsplot == "HIGH")] <- 18
	      plot(restab[, "N_sites"], restab[, "t_Entropy"], col = colsplot, pch = as.numeric(pchsplot), main = "Multilocus saturation test\nEntropy statistic", xlab = "N sites", ylab = "t (entropy statistic)")
	      legend("topleft", pch = c(NA, 16:18, NA, 16, 16), legend = risklegend, col = risklegcol)
	}
	if("t_CI" %in% colnames(restab)){
	      if(input$dataTreatment != "codonpos"){
	      		colsplot <- rep("black", nrow(restab))
			colsplot[which(restab[, "Risk_CI"] == "HIGH")] <- "red"
			colsplot[which(restab[, "Risk_CI"] == "LOW")] <- "blue"
              }
              pchsplot <- restab[, "Risk_CI"]
              pchsplot[which(pchsplot == "LOW")] <- 16
              pchsplot[which(pchsplot == "MEDIUM")] <- 17
              pchsplot[which(pchsplot == "HIGH")] <- 18
              plot(restab[, "N_sites"], restab[, "t_CI"], col = colsplot, pch = as.numeric(pchsplot), main = "Multilocus saturation test\nConsistency Index", xlab = "N sites", ylab = "t (consistency index)")
              legend("topleft", pch = c(NA, 16:18, NA, 16, 16), legend = risklegend, col = risklegcol)
        }
	if("t_Compression" %in% colnames(restab)){
	      if(input$dataTreatment != "codonpos"){
	      		colsplot <- rep("black", nrow(restab))
              		colsplot[which(restab[, "Risk_Compression"] == "HIGH")] <- "red"
            		colsplot[which(restab[, "Risk_Compression"] == "LOW")] <- "blue"
              }
              pchsplot <- restab[,"Risk_Compression"]
              pchsplot[which(pchsplot == "LOW")] <- 16
              pchsplot[which(pchsplot == "MEDIUM")] <- 17
              pchsplot[which(pchsplot == "HIGH")] <- 18
              plot(restab[, "N_sites"], restab[, "t_Compression"], col = colsplot, pch = as.numeric(pchsplot), main = "Multilocus saturation test\nCompression statistic", xlab = "N sites", ylab = "t (compression statistic)")
              legend("topleft", pch = c(NA, 16:18, NA, 16, 16), legend = risklegend, col = risklegcol)
        }
	dev.off()
}

setwd(initial.dir)

print("Assessment completed successfully.")