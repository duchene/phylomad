library(phangorn)

machine <- Sys.info()[["sysname"]]

ngene <- nrow(input$dataPath)

print(paste("The number of loci in this assessment was", ngene))


if(machine == "Darwin"){
	   iqtreePath <- paste0(getwd(), "/otherScripts/iqtree")
} else if(machine == "Linux"){
           iqtreePath <- paste0(getwd(), "/otherScripts/iqtree")
}

source("testBackbones/Substitutions/test.likelihood.multigene.R", local = T)