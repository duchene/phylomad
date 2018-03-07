source("otherScripts/R/clean.gene.R")
source("otherScripts/R/test.saturation.R")


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

geneResults <- test.saturation(loci = as.character(input$dataPath[, 4]), format = input$dataFormat, para = parallelise, ncore = input$Ncores, clean = input$dataTreatment, plotdat = plotdat)

#### Output missing

locinames <- as.character(input$dataPath[, 1])

if("tsat" %in% whatToOutput){
	
	#saturation.result <- 
	#write.csv(saturation.result, file = "saturation.test.phylomad")
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
	
}

setwd(initial.dir)

print("Assessment completed successfully.")