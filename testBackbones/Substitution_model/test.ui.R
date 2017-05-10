tabsetPanel(
	tabPanel("Data",
                       fluidRow(
		       column(1),
                       column(10,
		       h4("Select the data for which the model will be assessed"),
		       br(),
                       fileInput("dataPath", label = h5("Multiple data alignments in the same folder can also be selected"), multiple = T),
		       radioButtons("dataFormat", label = h4("Select the format of your data"),
                       choices = list("Phylip" = "phylip", "FASTA" = "fasta", "NEXUS" = "nexus"), selected = "phylip"),
		       fileInput("treesPath", label = h4("Select any input tree(s) (or leave blank)"), multiple = T),
                       radioButtons("treesFormat", label = h5("Select the format of your tree(s)"),
                       choices = list("No input tree" = "none", "NEWICK" = "newick", "NEXUS" = "nexus"), selected = "none")
		       ),
		       column(1))
        ),
	tabPanel("Model",
		       fluidRow(
		       column(1),
		       column(10,
		       selectInput("model", label = h4("Select the substitution model that will be assessed"),
		       choices = list("GTR+G" = "GTR+G", "GTR" = "GTR", "HKY+G" = "HKY+G", "HKY" = "HKY", "JC+G" = "JC+G", "JC" = "JC", "Automatic BIC model selection" = "autoModel"))
		       ),
		       column(1))
	),
	tabPanel("Test statistics",
		       fluidRow(
		       column(1),
		       column(10,
		       checkboxGroupInput("testStats", label = h4("Select test statistics for assessment"),
		       choices = list("Chi-Squared statistic" = "chisq", "Multinomial statistic" = "multlik", "Delta statistic" = "delta", "Biochemical diversity" = "biochemdiv", "Consistency index" = "consind", "Branch support" = "brsup", "Branch support 95% interval" = "CIbrsup", "Tree length" = "trlen", "Squared Mahalanobis distance" = "maha"), selected = c("chisq", "multlik", "delta", "biochemdiv", "consind", "brsup", "CIbrsup", "trlen", "maha")),
		       p("Note that the Mahalanobis statistic should only be selected if more than one other statistic is selected.")
		       ),
		       column(1))
	),
	tabPanel("Output",
		       fluidRow(
                       column(1),
                       column(10,
		       checkboxGroupInput("whatToOutput", label = h4("Select the output desired"),
                       choices = list("Metrics of adequacy" = "pvals", "Estimated tree for empirical data" = "phyloempres", "Simulated data" = "simdat", "Estimated trees for simulated data" = "phylosimres", "Test plots" = "testPlots"), selected = c("pvals", "testPlots")),
		       radioButtons("outputFormat", label = h4("Select the format of the output data"),
                       choices = list("Phylip" = "phylip", "FASTA" = "fasta", "NEXUS" = "nexus"), selected = "phylip"),
		       textInput("outputFolder", label = h4("Optional: Type the path of folder place output. Default is outputFolder in the main sofware folder"), value = paste0(getwd(), "/outputFolder/"))
		       ),
		       column(1))
	),
	tabPanel("Other options and START",
		       fluidRow(
                       column(1),
                       column(10,
		       numericInput("Nsims", label = h4("Select the number of simulations to be made"),
                       value = 100),
		       numericInput("Ncores", label = h4("Select the number of computer cores to be used for assessment"),
                       value = 1),
		       br(),
		       radioButtons("framework", label = h4("Select the statistical framework to use for assessment"),
                       choices = list("Likelihood" = "likelihood", "Bayesian" = "bayesian"), selected = "likelihood"),
		       br(),
		       actionButton("startAnalysis", label = "START ASSESSMENT"),
		       br(),
		       br(),
		       h5("Note that previous analyses in the output folder will be overwritten"),
		       h5("Progress is shown in the R console window"),
		       h5("The steps of analysis are shown below after completion"),
		       br(),
		       uiOutput("analysisConfirmation"),
		       br(),
		       uiOutput("analysisProgress")),
		       column(1))
	)
	
)