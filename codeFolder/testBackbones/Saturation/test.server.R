library(phangorn)

machine <- Sys.info()[["sysname"]]

ngene <- nrow(input$dataPath)

print(paste("The number of loci in this assessment was", ngene))

source("testBackbones/Saturation/test.saturation.multigene.R", local = T)