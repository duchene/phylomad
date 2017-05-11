library(phangorn)

machine <- Sys.info()[["sysname"]]

ngene <- nrow(input$dataPath)

print(paste("Number of loci and kind of machine have been identified as", ngene))

if(input$framework == "likelihood"){
	if(machine == "Darwin"){
		   phymlPath <- paste0(getwd(), "/otherScripts/PhyML-3.1/PhyML-3.1_macOS-MountainLion")
	} else if(machine == "Windows"){
		   phymlPath <- paste0(getwd(), "/otherScripts/PhyML-3.1/PhyML-3.1_win32.exe")
	} else if(machine == "Linux"){
	           phymlPath <- paste0(getwd(), "/otherScripts/PhyML-3.1/PhyML-3.1_linux64")
	}

	print(paste("Path for PhyML has been identified as", input$framework, "and machine as", machine))
	
	source("testBackbones/Substitution_model/test.likelihood.multigene.R", local = T)

} else if(input$framework == "bayesian"){
        if(machine == "Darwin"){
	           mbPath <-	paste0()
	} else if(machine == "Windows"){
		   mbPath <-	paste0()
	} else if(machine == "Linux"){
	       	   mbPath <-	paste0()
	}

	print("Path for PhyML has been identified")
	
	if(ngene == 1){
	      	source("testBackbones/Substitution_model/test.bayesian.gene.R", local = T)
	} else if(ngene > 1){
	        source("testBackbones/Substitution_model/test.bayesian.multigene.R", local = T)
	}
}