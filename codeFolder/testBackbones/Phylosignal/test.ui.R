tabsetPanel(
	tabPanel("Data",
                       fluidRow(
		       column(1),
                       column(10,
		       h4("Select the type of signal assessment"),
		       radioButtons("dataFormat", label = "",
                       choices = list("Locus (signal in sites)" = "locus", "Genome (signal in gene-trees)" = "genome"), selected = "locus"),
		       br(),
		       h4("Select the data set(s)"),
                       fileInput("dataPath", label = h5("Data sets must come from the same folder and be one or several alignments or trees files"), multiple = T),
		       br(),
		       h4("For alignments, select the format of your data"),
		       radioButtons("dataFormat", label = "",
                       choices = list("NEXUS" = "nexus", "NEWICK" = "newick", "Phylip" = "phylip", "FASTA" = "fasta"), selected = "nexus")
		       ),
		       column(1))
        ),
	tabPanel("Test statistics",
		       fluidRow(
		       column(1),
		       column(10,
		       tags$span(h4("Select test statistics for assessment"), tipify(bsButton("pB2", "?", style = "inverse", size = "extra-small"), placement = "right", "Test statistics are metrics used to compare empirical data with those generated under the model. They form the core of model assessment.")),
		       checkboxGroupInput("testStats", label = "",
		       choices = list("Dist to net-like point  " = "dnet", "Dist to tree-like point " = "dtree", "Entropy (ICA)  " = "entrop", "Internode cert (IC)  " = "icert", "Branch support  " = "brsup", "Binom P  " = "binp"), selected = c("dnet", "dtree", "entrop", "icert", "brsup", "binp")),
		       makeCheckboxTooltip(checkboxValue = "dnet", buttonLabel = "?", Tooltip = "Distance to closest network-like point of attraction in a ternary plot of the signals for each of the tree topologies proposed by a bipartition."),
		       makeCheckboxTooltip(checkboxValue = "dtree", buttonLabel = "?", Tooltip = "Distance to closest tree-like point of attraction in a ternary plot of the signals for each of the three topologies proposed by a bipartition."),
		       makeCheckboxTooltip(checkboxValue = "entrop", buttonLabel = "?", Tooltip = "Entropy in the signals of the topologies proposed by a bipartition."),
		       makeCheckboxTooltip(checkboxValue = "icert", buttonLabel = "?", Tooltip = "Entropy in the signals of the topologies proposed by a bipartition, using only the two highest values."),
		       makeCheckboxTooltip(checkboxValue = "brsup", buttonLabel = "?", Tooltip = "The mean node support across the inferred tree."),
		       makeCheckboxTooltip(checkboxValue = "binp", buttonLabel = "?", Tooltip = "P-value of binomial test of the superiority in the signal of the highest concordance factor in the data.")
		       ),
		       column(1))
	),
	tabPanel("Output",
		       fluidRow(
                       column(1),
                       column(10,
		       h4("Select the output desired"),
		       checkboxGroupInput("whatToOutput", label = "",
                       choices = list("Mean and per-branch P-values" = "pvals", "Test plots" = "testPlots", "Tree estimated from empirical data" = "phyloempres", "Simulated data" = "simdat", "All bipartition probabilities" = "allbipartp"), selected = c("pvals", "testPlots")),
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
		       checkboxInput("overwrite", label = h5("Overwrite previous analyses"), value = FALSE),
		       br(),
		       actionButton("startAnalysis", label = "START ASSESSMENT")
		       ),
		       column(1))
	)
	
)
