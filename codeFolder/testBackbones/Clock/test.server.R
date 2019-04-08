library(phangorn)
library(apTreeshape)

machine <- Sys.info()[["sysname"]]

ngene <- nrow(input$dataPath)

print(paste("The number of loci in this assessment was", ngene))

if(machine == "Darwin"){
	iqtreePath <- paste0(getwd(), "/otherScripts/iqtree")
} else if(machine == "Windows"){
	iqtreePath <- paste0(getwd(), "/otherScripts/iqtree")
} else if(machine == "Linux"){
	iqtreePath <- paste0(getwd(), "/otherScripts/iqtree")
}

source("testBackbones/Clock/test.clock.multigene.R", local = T)
