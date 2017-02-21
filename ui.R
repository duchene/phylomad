library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(

  # Application title
  titlePanel("PhyAdequacy"),

  sidebarLayout(
    sidebarPanel(
	h4("Welcome to PhyAdequacy"),
	p("This software assesses the adequacy of models used in phylogenetics.")

    ),
    mainPanel(
	tabsetPanel(
		tabPanel(
			fluidRow("Test type",
  	        		selectInput("modeltotest", label = h3("What model would you like to test, and in what statistical framework?"),
        			choices = as.list(dir("tests")), selected = 1)
			)
		),
		uiOutput("modeltestcode")
    )
  )
)))