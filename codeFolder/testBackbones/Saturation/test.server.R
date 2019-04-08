library(phangorn)

machine <- Sys.info()[["sysname"]]

if(machine == "Darwin"){
                   iqtreePath <- paste0(getwd(), "/otherScripts/iqtree")
} else if(machine == "Windows"){
                   iqtreePath <- paste0(getwd(), "/otherScripts/iqtree")
} else if(machine == "Linux"){
                   iqtreePath <- paste0(getwd(), "/otherScripts/iqtree")
}

ngene <- nrow(input$dataPath)

print(paste("The number of loci in this assessment was", ngene))

source("testBackbones/Saturation/test.saturation.multigene.R", local = T)