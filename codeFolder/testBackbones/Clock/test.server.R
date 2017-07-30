library(phangorn)
library(apTreeshape)

machine <- Sys.info()[["sysname"]]

ngene <- nrow(input$dataPath)

print(paste("The analysis contained", ngene, "loci"))

if(machine == "Darwin"){
	phymlPath <- paste0(getwd(), "/otherScripts/PhyML-3.1/PhyML-3.1_macOS-MountainLion")
} else if(machine == "Windows"){
	phymlPath <- paste0(getwd(), "/otherScripts/PhyML-3.1/PhyML-3.1_win32.exe")
} else if(machine == "Linux"){
	phymlPath <- paste0(getwd(), "/otherScripts/PhyML-3.1/PhyML-3.1_linux64")
}

print(paste("The machine type was identified as", machine))

source("testBackbones/Clock/test.clock.multigene.R", local = T)
