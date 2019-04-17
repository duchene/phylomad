library(phangorn)

machine <- Sys.info()[["sysname"]]

ngene <- nrow(input$dataPath)

print(paste("The number of loci in this assessment was", ngene))

iqtreePath <- paste0(getwd(), "/otherScripts/iqtree")
astralPath <- paste0(getwd(), "/otherScripts/Astral/astral.5.6.3.jar")

source("testBackbones/Phylosignal/test.likelihood.multigene.R", local = T)