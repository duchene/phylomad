library(phangorn)

machine <- Sys.info()[["sysname"]]

ngene <- nrow(input$dataPath)

dat <- readLines(as.character(input$dataPath$datapath[1]))

if(input$framework == "likelihood"){
	writeLines(dat, con = paste0(input$outputFolder, "pear.R"))
	#source("test.likelihood.R")
} else if(input$framework == "bayesian"){
        save(dat, file = paste0(input$outputFolder, "banana.R"))
	#source("test.bayesian.R")
}