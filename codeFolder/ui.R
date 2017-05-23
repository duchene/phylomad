library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(

  # Application title
  titlePanel(
  #h1("PhyloMAd", align = "center"),
  img(src = "phylomad.temp.png", height = 125, width = 900), #, style = "display: block; margin-left: auto; margin-right: auto;"),
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
    br(),
    h5("PhyloMAd is distributed under a yet unknown licence.")
    ),
    mainPanel(
    uiOutput("modeltestcode")
  )
)))