library(phangorn)

machine <- Sys.info()[["sysname"]]

if(machine == "Darwin"){
                   phymlPath <- paste0(getwd(), "/otherScripts/PhyML-3.1/PhyML-3.1_macOS-MountainLion")
} else if(machine == "Windows"){
                   phymlPath <- paste0(getwd(), "/otherScripts/PhyML-3.1/PhyML-3.1_win32.exe")
} else if(machine == "Linux"){
                   phymlPath <- paste0(getwd(), "/otherScripts/PhyML-3.1/PhyML-3.1_linux64")
}


ngene <- nrow(input$dataPath)

print(paste("The number of loci in this assessment was", ngene))

source("testBackbones/Saturation/test.saturation.multigene.R", local = T)