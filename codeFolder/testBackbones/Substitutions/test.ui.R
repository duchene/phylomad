tabsetPanel(
	tabPanel("Data",
                       fluidRow(
		       column(1),
                       column(10,
		       h4("Select the nucleotide or amino acid alignment(s) for which the substitution model will be assessed"),
                       fileInput("dataPath", label = h5("Alignments selected must come from the same folder"), multiple = T),
		       radioButtons("dataFormat", label = h4("Select the format of your data"),
                       choices = list("NEXUS" = "nexus", "Phylip" = "phylip", "FASTA" = "fasta"), selected = "nexus"),
		       fileInput("treesPath", label = h4("Select any input tree(s) (or leave blank)"), multiple = T),
                       radioButtons("treesFormat", label = h5("Select the format of your tree(s)"),
                       choices = list("No input tree" = "none", "NEXUS" = "nexus", "NEWICK" = "newick"), selected = "none")
		       ),
		       column(1))
        ),
	tabPanel("Model",
		       fluidRow(
		       column(1),
		       column(10,
		       selectInput("model", label = h4("Select the substitution model that will be assessed"),
		       choices = list("GTR" = "GTR", "HKY" = "HKY", "JC" = "JC", "Automatic BIC model selection" = "autoModel", "JTT" = "JTT", "LG" = "LG", "WAG" = "WAG", "Dayhoff" = "Dayhoff")),
		       radioButtons("RASmodel", label = h4("Select the model of rates across sites"),
                       choices = list("Gamma-distributed" = "+G", "Equal across sites" = ""), selected = "+G")
		       ),
		       column(1))
	),
	tabPanel("Test statistics",
		       fluidRow(
		       column(1),
		       column(10,
		       checkboxGroupInput("testStats", label = h4("Select test statistics for assessment"),
		       choices = list("Chi-squared statistic" = "chisq", "Multinomial statistic" = "multlik", "Delta statistic" = "delta", "Biochemical diversity" = "biochemdiv", "Consistency index" = "consind", "Branch support" = "brsup", "Branch support 95% interval" = "CIbrsup", "Tree length" = "trlen", "Squared Mahalanobis distance" = "maha"), selected = c("chisq", "multlik", "delta", "biochemdiv", "consind", "brsup", "CIbrsup", "trlen", "maha")),
		       p("Only select the Mahalanobis distance if multiple other statistics are also selected.")
		       ),
		       column(1))
	),
	tabPanel("Output",
		       fluidRow(
                       column(1),
                       column(10,
		       checkboxGroupInput("whatToOutput", label = h4("Select the output desired"),
                       choices = list("Only overall summary file (overrides other options)" = "simple", "Metrics of adequacy for individual loci" = "pvals", "Tree estimated from empirical data" = "phyloempres", "Simulated data" = "simdat", "Trees estimated from simulated data" = "phylosimres", "Test plots" = "testPlots"), selected = c("pvals", "testPlots")),
		       radioButtons("outputFormat", label = h4("Select the format of the output data"),
                       choices = list("NEXUS" = "nexus", "Phylip" = "phylip", "FASTA" = "fasta"), selected = "nexus"),
		       textInput("outputFolder", label = h4("Optional: Type the path of the output folder. The default is outputFolder in the main sofware folder"), value = paste0(getwd(), "/../outputFolder/"))
		       ),
		       column(1))
	),
	tabPanel("Other options and START",
		       fluidRow(
                       column(1),
                       column(10,
		       numericInput("Nsims", label = h4("Select the number of simulations to be made"),
                       value = 100),
		       h4("Select the number of computer cores to be used"),
		       numericInput("Ncores", label = h5("Multi-core assessments can only be aborted at the completion of a locus assessment"),
                       value = 1),
		       br(),
		       radioButtons("framework", label = h4("Select the statistical framework to use for assessment"),
                       choices = list("Likelihood" = "likelihood", "Bayesian (temporarily unavailable)" = "bayesian"), selected = "likelihood"),
		       br(),
		       actionButton("startAnalysis", label = "START ASSESSMENT"),
		       br(),
		       br(),
		       h5("OUTPUT NOTES"),
		       h5("Previous results with the same paths will be overwritten"),
		       h5("Progress is shown in the shell window (mac) or in a log file saved to the main PhyloMAd folder (windows)"), 
		       h5("To avoid analyses running repeatedly, only press the start button once and be aware that the analyses might take a few minutes to print any progress"),
		       h5("Further notes, warnings, or errors will appear below after completion"),
		       br(),
		       uiOutput("analysisConfirmation"),
		       br(),
		       uiOutput("analysisProgress")),
		       column(1))
	)
	
)
