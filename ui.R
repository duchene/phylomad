library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(

  # Application title
  titlePanel(
  #h1("PhyloMAd", align = "center"),
  img(src = "phylomad.temp.png", height = 150, width = 300, style = "display: block; margin-left: auto; margin-right: auto;"),
  windowTitle = "PhyloMAd"
  ),
  sidebarLayout(
    sidebarPanel(
    h3("Software for assessment of phylogenetic model adequacy", align = "center"),
    h4("Start here", align = "center"),
    br(),
    radioButtons("modeltotest", label = h4("Select the model to be assessed"), 
    choices = rev(as.list(dir("testBackbones"))), selected = "Substitution_model"),
    br(),
    p("Note that several R packages other than shiny are needed to run model assessment."),
    p("Press the following button if you would like to install these packages. This might take a couple of minutes, and will complete with a new message below."),
    br(),
    actionButton("installPackages", label = "Install required packages"),
    br(),
    uiOutput("installationConfirmation")
    ),
    mainPanel(
    uiOutput("modeltestcode")
  )
)))