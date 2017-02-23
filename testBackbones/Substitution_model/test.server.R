library(phangorn)

machine <- Sys.info()[["sysname"]]

ngene <- nrow(input$dataPath)

if(input$framework == "likelihood"){
	if(machine == "Darwin"){
		   phymlPath <- paste0(getwd(), "/otherScripts/PhyML-3.1/PhyML-3.1_macOS-MountainLion")
	} else if(machine == "Windows"){
		   phymlPath <- paste0(getwd(), "/otherScripts/PhyML-3.1/PhyML-3.1_win32.exe")
	} else if(machine == "Linux"){
	           phymlPath <- paste0(getwd(), "/otherScripts/PhyML-3.1/PhyML-3.1_linux64")
	}
	if(ngene == 1){
		source("test.likelihood.gene.R")
	} else if(ngene > 1){
	       	source("test.likelihood.multigene.R")
	}
} else if(input$framework == "bayesian"){
        if(machine == "Darwin"){
	           mbPath <-	paste0()
	} else if(machine == "Windows"){
		   mbPath <-	paste0()
	} else if(machine == "Linux"){
	       	   mbPath <-	paste0()
	}
	if(ngene == 1){
	      	source("test.bayesian.gene.R")
	} else if(ngene > 1){
	        source("test.bayesian.multigene.R")
	}
}