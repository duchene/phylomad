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
	fluidRow(
		column(3),
		column(3,
      selectInput("Select your test", label = h3("Select your test"), 
        choices = list("Choice 1" = 1, "Choice 2" = 2,
                       "Choice 3" = 3), selected = 1)),
		column(3)


	)
    )
  )
))