library(phangorn)

machine <- Sys.info()[["sysname"]]

ngene <- nrow(input$dataPath)

print(paste("The number of loci in this assessment was", ngene))

if(input$framework == "likelihood"){
	if(machine == "Darwin"){
		   phymlPath <- paste0(getwd(), "/otherScripts/PhyML-3.1/PhyML-3.1_macOS-MountainLion")
	} else if(machine == "Windows"){
		   phymlPath <- paste0(getwd(), "/otherScripts/PhyML-3.1/PhyML-3.1_win32.exe")
	} else if(machine == "Linux"){
	           phymlPath <- paste0(getwd(), "/otherScripts/PhyML-3.1/PhyML-3.1_linux64")
	}

	print(paste("Analysis was performed in the", input$framework, "statistical framework and the machine type is", if(machine == "Darwin") "mac" else machine))
	
	source("testBackbones/Substitutions/test.likelihood.multigene.R", local = T)

} else if(input$framework == "bayesian"){
        if(machine == "Darwin"){
	           mbPath <-	paste0()
	} else if(machine == "Windows"){
		   mbPath <-	paste0()
	} else if(machine == "Linux"){
	       	   mbPath <-	paste0()
	}

	print(paste("Analysis was performed in the", input$framework, "statistical framework and the machine type is", machine))
	
	if(ngene == 1){
	      	source("testBackbones/Substitutions/test.bayesian.gene.R", local = T)
	} else if(ngene > 1){
	        source("testBackbones/Substitutions/test.bayesian.multigene.R", local = T)
	}
}