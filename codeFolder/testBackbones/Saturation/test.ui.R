tabsetPanel(
	tabPanel("Data",
                       fluidRow(
		       column(1),
                       column(10,
		       h4("Select the nucleotide or amino acid alignment(s) for which the substitution model will be assessed"),
                       fileInput("dataPath", label = h5("Alignments selected must come from the same folder"), multiple = T),
		       br(),
		       h4("Select the method of data handling"),
		       radioButtons("dataTreatment", label = "", choices = list("Remove sites with missing data  " = "cleandata", "Test complete data set  " = "complete", "1st+2nd versus 3rd codon positions  " = "codonpos")),
		       makeCheckboxTooltip(checkboxValue = "cleandata", buttonLabel = "?", Tooltip = "This will lead to the most accurate test of saturation, but might be a poor reflection of the quality of phylogenetic inferences if a large amount of data are missing."),
                       makeCheckboxTooltip(checkboxValue = "complete", buttonLabel = "?", Tooltip = "Loci will be kept whole for assessment. This can lead to a misleading saturation test if there is a large amount of missing data in the alignment."),
                       makeCheckboxTooltip(checkboxValue = "codonpos", buttonLabel = "?", Tooltip = "This performs independent tests of saturation for the two types of sites.")
		       ),
		       column(1))
        ),
	tabPanel("Output",
		       fluidRow(
		       column(1),
		       column(10,
		       h4("Select saturation statistics"),
                       checkboxGroupInput("saturationStats", label = "",
                       choices = list("Entropy" = "enth", "Consistency Index" = "cith", "Compression statistic" = "comth"), selected = c("enth", "cith", "comth")),
                       br(),
		       h4("Select the output desired"),
                       checkboxGroupInput("whatToOutput", label = "",
                       choices = list("Test of saturation for loci" = "tsat", "Saturation plots" = "satPlots", "Multi-locus saturation plots" = "multiSatPlots"), selected = c("tsat", "satPlots")),
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
