library(phangorn)
library(apTreeshape)

machine <- Sys.info()[["sysname"]]

ngene <- nrow(input$dataPath)

print(paste("Number of loci and kind of machine have been identified as", ngene))

if(machine == "Darwin"){
	phymlPath <- paste0(getwd(), "/otherScripts/PhyML-3.1/PhyML-3.1_macOS-MountainLion")
} else if(machine == "Windows"){
	phymlPath <- paste0(getwd(), "/otherScripts/PhyML-3.1/PhyML-3.1_win32.exe")
} else if(machine == "Linux"){
	phymlPath <- paste0(getwd(), "/otherScripts/PhyML-3.1/PhyML-3.1_linux64")
}

print(paste("Path for PhyML has been identified as", input$framework, "and machine as", machine))

source("testBackbones/Clock_model/test.clock.multigene.R", local = T)
