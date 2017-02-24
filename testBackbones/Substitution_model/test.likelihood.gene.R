source("otherScripts/R/run.gene.R")
source("otherScripts/R/get.test.statistics.R")
source("otherScripts/R/runPhyML.R")
source("otherScripts/R/clean.gene.R")
source("otherScripts/R/get.model.R")
if("chisq" %in% unlist(input$testStats)) source("testStatistics/get.chisqstat.R")
if("biochemdiv" %in% unlist(input$testStats)) source("testStatistics/get.biodivstat.R")

print("Functions required have been loaded")

if(input$model == "autoModel") modeltested <- get.model(as.character(input$dataPath[1, 4])) else modeltested <- input$model

print("Model has been identified")

geneDNAbin <- clean.gene(sdata = as.character(input$dataPath[1, 4]), format = input$dataFormat)

print("Gene has been cleaned")

if(input$Ncores > 1) parallelise <- T else parallelise <- F

setwd(input$outputFolder)

print("Output folder has been identified")

geneResults <- run.gene(sdata = geneDNAbin, format = "DNAbin", model = modeltested, phymlPath = phymlPath, Nsims = input$Nsims, para = parallelise, ncore = input$Ncores, testStats = unlist(input$testStats))

print("Gene results have been processed")

if("pvals" %in% unlist(input$whatToOutput)){
	out <- rbind(unlist(geneResults[grep("[.]tailp", names(geneResults))]), unlist(geneResults[grep("emp[.]", names(geneResults))]), unlist(geneResults[grep("[.]sdpd", names(geneResults))]))
	write.csv(out, file = "output.pvals.PhyloMAd.csv")
}

if("simdat" %in% unlist(input$whatToOutput)){
        
}

if("testPlots" %in% unlist(input$whatToOutput)){
	
}
