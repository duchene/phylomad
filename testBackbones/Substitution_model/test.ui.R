tabsetPanel(
	tabPanel("Data",
                       fluidRow(
                       column(3),
                       column(6,
                       fileInput("dataPath", label = h4("Type the path to the alignment or folder with alignments for which the model will be assessed", align = "center")),
		       radioButtons("locusOrGenome", label = h4("Select whether your data are a single or multiple alignments", align = "center"),
                       choices = list("Single" = "single", "Multiple" = "multiple"), selected = "single")
		       ),
                       column(3))
        ),
	tabPanel("Model",
		       fluidRow(
		       column(3),
		       column(6,
		       selectInput("model", label = h4("Select the substitution model that will be assessed", align = "center"),
		       choices = list("GTR+G" = "GTR+G", "GTR" = "GTR", "HKY+G" = "HKY+G", "HKY" = "HKY", "JC+G" = "JC+G", "JC" = "JC")),
		       numericInput("Nsims", label = h4("Select the number of simulations to be made", align = "center"),
                       value = 100)
		       ),
		       column(3))
	),
	tabPanel("Test statistics",
		       fluidRow(
		       column(3),
		       column(6,
		       checkboxGroupInput("testStats", label = h4("Select test statistics for assessment", align = "center"),
		       choices = list("Chi-Squared statistic" = 1, "Multinomial statistic" = 2, "Delta" = 3, "Biochemical diversity" = 4, "Consistency index" = 5, "Branch support" = 6, "Tree length" = 7, "Mahalanobis" = 8), selected = 1
		       ),
		       column(3)))
	),
	tabPanel("Output",
		       fluidRow(
                       column(3),
                       column(6,
		       checkboxGroupInput("whatToOutput", label = h4("Select the output desired", align = "center"),
                       choices = list("P-values" = 1, "Simulated data" = 2, "Test plots" = 3), selected = 1),
		       textInput("outputFolder", label = h4("Full path of folder where output should be saved", align = "center"), 
	               value = getwd()),
		       column(3)))
	),
	tabPanel("Other options and START",
		       fluidRow(
                       column(3),
                       column(6,
		       radioButtons("machine", label = h4("Select the kind of machine that you are using", align = "center"),
        	       choices = list("Linux" = "linux", "Mac" = "mac", "Windows" = "windows"), selected = "mac"),
		       numericInput("Ncores", label = h4("Select the number of computer cores to be used for assessment", align = "center"),
                       value = 1),
		       br(),
		       radioButtons("framework", label = h4("Choose the statistical framework to use for assessment", align = "center"),
                       choices = list("Likelihood" = "likelihood", "Bayesian" = "bayesian"), selected = "likelihood"),
		       br(),
		       actionButton("startAnalysis", label = "START ASSESSMENT"),
		       br(),
		       p("After the analysis has started, there will be two messages below to indicate progress, and a final message to indicate completion. You can also look at your R console to observe some of the progress."),
		       br(),
		       uiOutput("progressMessage1"),
		       uiOutput("progressMessage2"),
		       uiOutput("analysisConfirmation")),
		       column(3))
	)
	
)
