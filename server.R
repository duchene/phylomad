library(shiny)

options(shiny.maxRequestSize = 100*1024^2)

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {

  output$files = renderDataTable({
    files = list.files(readDirectoryInput(session, "directory"), full.names = T)
    data.frame(name = basename(files), file.info(files))
  })

  output$modeltestcode <- renderUI({
  	source(paste0("testBackbones/", input$modeltotest, "/test.ui.R"), local = T)  
  })
  
  observeEvent(input$installPackages, {
  	source(paste0("testBackbones/", input$modeltotest, "/install.test.packages.R"), local = T)
        output$installationConfirmation <- renderText("The packages are now installed")
  })

  observeEvent(input$startAnalysis, { 
  	output$analysisProgress <- renderPrint(source(paste0("testBackbones/", input$modeltotest, "/test.server.R"), local = T))
	output$analysisConfirmation <- renderText("The analysis has been completed")
  })

  session$onSessionEnded(function() { 
     stopApp()
     #q("no") 
  })

  })