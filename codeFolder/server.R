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

  observeEvent(input$startAnalysis, {
  		       showModal(modalDialog(
				title = "Are you sure?",
				"This will start model assessment. Progress will be shown in the terminal window (mac), or in the log file saved in the main PhyloMAd folder (windows).

				Previous results with the same folder and file names will be overwritten.",
				easyClose = F,
				footer = tagList(
				       modalButton("Cancel"),
          			       actionButton("START", "START")
        			)
		      ))
  })

  observeEvent(input$START, {
  		       showModal(modalDialog(
				title = "The steps undergone, warnings, and errors, will appear below on completion.\n\nPressing 'Close' will not stop the analysis (only closing PhyloMAd will).\n\nIf you close this window, avoid pressing START ANALYSIS again, since this will run several analyses repeatedly.",
				renderPrint(source(paste0("testBackbones/", input$modeltotest, "/test.server.R"), local = T)),
				easyClose = F,
				footer = modalButton("Close")
		       ))
  })

  session$onSessionEnded(function() { 
     stopApp()
     q("no") 
  })

  })