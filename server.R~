library(shiny)

options(shiny.maxRequestSize = 100*1024^2)

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {

# The following is to read the new directoryInput function
  observeEvent(
  ignoreNULL = TRUE,
  eventExpr = {
    input$directory
  },
  handlerExpr = {
    if (input$directory > 0) {
      path = choose.dir(default = readDirectoryInput(session, "directory"))
      updateDirectoryInput(session, "directory", value = path)
    }
  }
  )

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
  	source(paste0("testBackbones/", input$modeltotest, "/test.server.R"), local = T)
	output$analysisConfirmation <- renderText("The analysis has been completed")
  })

  })