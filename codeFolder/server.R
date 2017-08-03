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
				h5("This will start model assessment."),
				h5("Progress will be shown in the terminal window (mac), or in the log file saved in the main PhyloMAd folder (windows)."),
				easyClose = F,
				footer = tagList(
				       modalButton("Cancel"),
          			       actionButton("START", "START")
        			)
		      ))
  })

  observeEvent(input$START, {
  		       showModal(modalDialog(
				title = "Assessment steps appear below on completion.",
				h5("'Close' will not stop the analysis (closing PhyloMAd will)."),
				h5("Avoid pressing START ANALYSIS again, since analyses will be run repeatedly."),
				h5("In the event of an error, PhyloMAd should be closed and reopened."),
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