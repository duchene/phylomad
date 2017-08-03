tabsetPanel(
	tabPanel("Data",
                       fluidRow(
		       column(1),
                       column(10,
		       h4("Select the nucleotide or amino acid alignment(s) for which the substitution model will be assessed"),
                       fileInput("dataPath", label = h5("Alignments selected must come from the same folder"), multiple = T),
		       br(),
		       h4("Select the format of your data"),
		       radioButtons("dataFormat", label = "",
                       choices = list("NEXUS" = "nexus", "Phylip" = "phylip", "FASTA" = "fasta"), selected = "nexus"),
		       br(),
		       h4("Select any input tree(s) (or leave blank)"),
		       fileInput("treesPath", label = "", multiple = T),
		       br(),
		       h4("Select the format of your tree(s)"),
                       radioButtons("treesFormat", label = "",
                       choices = list("No input tree" = "none", "NEXUS" = "nexus", "NEWICK" = "newick"), selected = "none")
		       ),
		       column(1))
        ),
	tabPanel("Model",
		       fluidRow(
		       column(1),
		       column(10,
		       h4("Select the substitution model that will be assessed"),
		       selectInput("model", label = "",
		       choices = list("Nucleotide - GTR" = "GTR", "Nucleotide - HKY" = "HKY", " Nucleotide - JC" = "JC", "Automatic nucleotide model selection (BIC)" = "autoModel", "Amino acid - JTT" = "JTT", "Amino acid - LG" = "LG", "Amino acid - WAG" = "WAG", "Amino acid - Dayhoff" = "Dayhoff")),
		       br(),
		       h4("Select the model of rates across sites"),
		       radioButtons("RASmodel", label = "",
                       choices = list("Gamma-distributed" = "+G", "Equal across sites" = ""), selected = "+G")
		       ),
		       column(1))
	),
	tabPanel("Test statistics",
		       fluidRow(
		       column(1),
		       column(10,
		       h4("Select test statistics for assessment"),
		       checkboxGroupInput("testStats", tags$span("", tipify(bsButton("pB2", "?", style = "inverse", size = "extra-small"), "Hover over each test statistic for further details.")),
		       choices = list("Chi-squared statistic" = "chisq", "Multinomial statistic" = "multlik", "Delta statistic" = "delta", "Biochemical diversity" = "biochemdiv", "Consistency index" = "consind", "Branch support" = "brsup", "Branch support 95% interval" = "CIbrsup", "Tree length" = "trlen", "Squared Mahalanobis distance" = "maha"), selected = c("chisq", "multlik", "delta", "biochemdiv", "consind", "brsup", "CIbrsup", "trlen", "maha")),
		       makeCheckboxTooltip(checkboxValue = "chisq", buttonLabel = "?", Tooltip = "To be added"),
		       makeCheckboxTooltip(checkboxValue = "multlik", buttonLabel = "?", Tooltip = "To be added"),
		       makeCheckboxTooltip(checkboxValue = "delta", buttonLabel = "?", Tooltip = "To be added"),
		       makeCheckboxTooltip(checkboxValue = "biochemdiv", buttonLabel = "?", Tooltip = "To be added"),
		       makeCheckboxTooltip(checkboxValue = "consind", buttonLabel = "?", Tooltip = "To be added"),
		       makeCheckboxTooltip(checkboxValue = "brsup", buttonLabel = "?", Tooltip = "To be added"),
		       makeCheckboxTooltip(checkboxValue = "CIbrsup", buttonLabel = "?", Tooltip = "To be added"),
		       makeCheckboxTooltip(checkboxValue = "trlen", buttonLabel = "?", Tooltip = "To be added"),
		       makeCheckboxTooltip(checkboxValue = "maha", buttonLabel = "?", Tooltip = "To be added")
		       ),
		       column(1))
	),
	tabPanel("Output",
		       fluidRow(
                       column(1),
                       column(10,
		       h4("Select the output desired"),
		       checkboxGroupInput("whatToOutput", label = "",
                       choices = list("Metrics of adequacy for individual loci" = "pvals", "Test plots" = "testPlots", "Only overall summary file (overrides other options)" = "simple", "Tree estimated from empirical data" = "phyloempres", "Simulated data" = "simdat", "Trees estimated from simulated data" = "phylosimres"), selected = c("pvals", "testPlots")),
		       br(),
		       h4("Select the format of the output data"),
		       radioButtons("outputFormat", label = "",
                       choices = list("NEXUS" = "nexus", "Phylip" = "phylip", "FASTA" = "fasta"), selected = "nexus"),
		       br(),
		       h4("Optional: Type the path of the output folder. The default is outputFolder in the main sofware folder"),
		       textInput("outputFolder", label = "", value = paste0(getwd(), "/../outputFolder/"))
		       ),
		       column(1))
	),
	tabPanel("Other options and START",
		       fluidRow(
                       column(1),
                       column(10,
		       h4("Select the number of simulations to be made"),
		       numericInput("Nsims", label = "",
                       value = 100),
		       br(),
		       h4("Select the number of computer cores to be used"),
		       numericInput("Ncores", label = h5("Multi-core assessments can only be aborted at the completion of a locus assessment"),
                       value = 1),
		       br(),
		       h4("Select the statistical framework to use for assessment"),
		       radioButtons("framework", label = "",
                       choices = list("Likelihood" = "likelihood", "Bayesian (not available yet)" = "bayesian"), selected = "likelihood"),
		       br(),
		       checkboxInput("overwrite", label = h5("Overwrite previous analyses"), value = FALSE),
		       br(),
		       actionButton("startAnalysis", label = "START ASSESSMENT")
		       ),
		       column(1))
	)
	
)
