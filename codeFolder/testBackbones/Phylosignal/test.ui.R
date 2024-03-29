tabsetPanel(
	tabPanel("Data",
                       fluidRow(
		       column(1),
                       column(10,
		       br(),
		       h4("Select the type of assessment"),
		       radioButtons("testType", label = "",
		       # Option test gene-trees = "genome" is not currently here, but uses Astral and is meant to assess the assumption of even discordance factors under ILS across gene trees.
                       choices = list("Loci" = "locus", "Species tree" = "tree"), selected = "locus"),
		       makeCheckboxTooltip(checkboxValue = "locus", buttonLabel = "?", Tooltip = "Provides statstics on each of a set of loci, and does not require an input species tree."),
                       makeCheckboxTooltip(checkboxValue = "tree", buttonLabel = "?", Tooltip = "Only provides statistics to an input tree hypothesis, based on an input tree and loci."),
		       br(),
		       h4("Select the per-locus data"),
                       fileInput("dataPath", label = h5("Data sets must come from the same folder and be one or several alignments or trees files"), multiple = T),
		       checkboxInput("dataType", label = "Amino acid data", value = FALSE),
		       br(),
		       h4("Select tree to be tested"),
                       fileInput("sptreePath", label = h5("Only required if using the Test tree option"), multiple = T)
		       ),
		       column(1))
        ),
	tabPanel("Test statistics",
		       fluidRow(
		       column(1),
		       column(10,
		       tags$span(h4("Select test statistics for assessment"), tipify(bsButton("pB2", "?", style = "inverse", size = "extra-small"), placement = "right", "Test statistics are metrics used to compare empirical data with those generated under the model. They form the core of model assessment.")),
		       checkboxGroupInput("testStats", label = "",
		       choices = list("Dist to net-like point  " = "dnet", "Dist to tree-like point " = "dtree", "Entropy (ICA)  " = "entrop", "Internode cert (IC)  " = "icert", "Concordance factor  " = "CF", "Binom P  " = "binp", "D-statistic  " = "dstat", "KC statistic" = "kcstat"), selected = c("dnet", "dtree", "entrop", "icert", "CF", "binp", "dstat", "kcstat")),
		       makeCheckboxTooltip(checkboxValue = "dnet", buttonLabel = "?", Tooltip = "Distance to closest network-like point of attraction in a ternary plot of the signals for each of the tree topologies proposed by a quartet."),
		       makeCheckboxTooltip(checkboxValue = "dtree", buttonLabel = "?", Tooltip = "Distance to closest tree-like point of attraction in a ternary plot of the signals for each of the three topologies proposed by a quartet."),
		       makeCheckboxTooltip(checkboxValue = "entrop", buttonLabel = "?", Tooltip = "Entropy in the signals of the topologies proposed by a quartet."),
		       makeCheckboxTooltip(checkboxValue = "icert", buttonLabel = "?", Tooltip = "Entropy in the signals of the topologies proposed by a quartet, using only the two highest values."),
		       makeCheckboxTooltip(checkboxValue = "CF", buttonLabel = "?", Tooltip = "The proportional support of alignment sites or genes for the estimated quartet."),
		       makeCheckboxTooltip(checkboxValue = "binp", buttonLabel = "?", Tooltip = "P-value of binomial test of the superiority in the signal of the highest concordance factor in the data."),
		       makeCheckboxTooltip(checkboxValue = "dstat", buttonLabel = "?", Tooltip = "Invariants-based function of the two least represented trees in a quartet."),
                       makeCheckboxTooltip(checkboxValue = "kcstat", buttonLabel = "?", Tooltip = "Invariants-based function of the three quartet proportions.")
		       ),
		       column(1))
	),
	tabPanel("Output",
		       fluidRow(
                       column(1),
                       column(10,
		       h4("Select the output desired"),
		       checkboxGroupInput("whatToOutput", label = "",
                       choices = list("Tree with branch statistics" = "phyloempres", "P-values and bootstrap table" = "pvals", "Test plots" = "testPlots", "Simulated sequences" = "simdat", "Discordance factors" = "allqp"), selected = c("phyloempres", "pvals", "testPlots")),
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
                       numericInput("Ncores", label = h5("Multi-core assessments can only be aborted at the completion of a locus assessment. Only used for computing support metrics for an input species tree."), value = 1),
		       br(),
		       checkboxInput("overwrite", label = h5("Overwrite previous analyses"), value = FALSE),
		       br(),
		       actionButton("startAnalysis", label = "START ASSESSMENT")
		       ),
		       column(1))
	)
	
)
