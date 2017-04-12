tabsetPanel(
	tabPanel("Data",
                       fluidRow(
                       column(3),
                       column(6,
                       fileInput("dataPath", label = h4("Select the data for which the model will be assessed", align = "center"), multiple = T),
		       h5("You can select a single or multiple data alignments in the same folder."),
		       radioButtons("dataFormat", label = h4("Select the format of your data", align = "center"),
                       choices = list("Phylip" = "phylip", "FASTA" = "fasta", "NEXUS" = "nexus"), selected = "phylip")
                       ),
                       column(3))
        ),
	tabPanel("Model",
		       fluidRow(
		       column(3),
		       column(6,
		       selectInput("model", label = h4("Select the substitution model that will be assessed", align = "center"),
		       choices = list("GTR+G" = "GTR+G", "GTR" = "GTR", "HKY+G" = "HKY+G", "HKY" = "HKY", "JC+G" = "JC+G", "JC" = "JC", "Automatic BIC model selection" = "autoModel"))
		       ),
		       column(3))
	),
	tabPanel("Test statistics",
		       fluidRow(
		       column(3),
		       column(6,
		       checkboxGroupInput("testStats", label = h4("Select test statistics for assessment", align = "center"),
		       choices = list("Chi-Squared statistic" = "chisq", "Multinomial statistic" = "multlik", "Delta statistic" = "delta", "Biochemical diversity" = "biochemdiv", "Consistency index" = "consind", "Branch support" = "brsup", "Branch support 95% interval" = "CIbrsup", "Tree length" = "trlen", "Squared Mahalanobis distance" = "maha"), selected = c("chisq", "multlik", "delta", "biochemdiv", "consind", "brsup", "CIbrsup", "trlen", "maha")),
		       p("Note that the Mahalanobis statistic should only be selected if more than one other statistic is selected.")
		       ),
		       column(3))
	),
	tabPanel("Output",
		       fluidRow(
                       column(3),
                       column(6,
		       checkboxGroupInput("whatToOutput", label = h4("Select the output desired", align = "center"),
                       choices = list("Metrics of adequacy" = "pvals", "Phylogenetic analyses for empirical data" = "phyloempres", "Simulated data" = "simdat", "Phylogenetic analyses for simulated data" = "phylosimres", "Test plots" = "testPlots"), selected = c("pvals", "testPlots")),
		       radioButtons("outputFormat", label = h4("Select the format of the output data", align = "center"),
                       choices = list("Phylip" = "phylip", "FASTA" = "fasta", "NEXUS" = "nexus"), selected = "phylip"),
		       textInput("outputFolder", label = h4("Type in the path to folder in which to save the output. PhyloMAd offers a default outputFolder, in the main sofware folder"), value = paste0(getwd(), "/outputFolder/"))
		       ),
		       column(3))
	),
	tabPanel("Other options and START",
		       fluidRow(
                       column(3),
                       column(6,
		       numericInput("Nsims", label = h4("Select the number of simulations to be made", align = "center"),
                       value = 100),
		       numericInput("Ncores", label = h4("Select the number of computer cores to be used for assessment", align = "center"),
                       value = 1),
		       br(),
		       radioButtons("framework", label = h4("Select the statistical framework to use for assessment", align = "center"),
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
