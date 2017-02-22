library(shiny)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {

  output$modeltestcode <- renderUI({
  	source(paste0("testBackbones/", input$modeltotest, "/test.ui.R"), local = T)  
  })
  
  observeEvent(input$installPackages, {
  	source(paste0("testBackbones/", input$modeltotest, "/install.test.packages.R"), local = T)
        output$installationConfirmation <- renderText("The packages are now installed")
  })

  observeEvent(input$startAnalysis, { 
  	source(paste0("testBackbones/", input$modeltotest, "/test.server.R"), local = T)
	output$analysisConfirmation <- renderText("The analysis has been completed")
  })

  })