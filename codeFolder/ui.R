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
    h6("Copyright 2017 - PhyloMAd authors", align = "center"),
    h3("Software for Assessment of Phylogenetic Model Adequacy", align = "center"),
    br(),
    radioButtons("modeltotest", label = h4("Select the type of model to be assessed"), 
    choices = rev(as.list(dir("testBackbones"))), selected = "Substitutions")
    ),
    mainPanel(
    uiOutput("modeltestcode")
  )
)))