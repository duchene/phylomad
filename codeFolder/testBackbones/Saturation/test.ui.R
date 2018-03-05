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
		       checkboxInput("cleanOrNot", label = tags$span(h5("Remove sites with missing data from each alignment"), tipify(bsButton("pB1", "?", style = "inverse", size = "extra-small"), placement = "right", "Removing sites with missing data is important for accurate estimation of test statistics. Only uncheck this box if the alignments only have a negligible amount of missing data.")), value = TRUE)
		       ),
		       column(1))
        ),
	tabPanel("Output and START",
		       fluidRow(
                       column(1),
                       column(10,
		       h4("Select the output desired"),
		       checkboxGroupInput("whatToOutput", label = "",
                       choices = list("Test of saturation for loci" = "tsat", "Saturation plots" = "satPlots", "Multi-locus saturation plots" = "multiSatPlots"), selected = c("tsat", "satPlots")),
		       br(),
		       h4("Optional: Type the path of the output folder. The default is outputFolder in the main sofware folder"),
		       textInput("outputFolder", label = "", value = paste0(getwd(), "/../outputFolder/")),
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
