get.test.statistics <- function(sdata, format = "phylip", aadata = F, geneName = "empirical", iqtreePath, model = "GTR+G", stats = c("chisq", "multlik", "delta", "biochemdiv", "consind", "brsup", "trlen", "maha"), tree = NULL, getTreeForced = F, ncore = 1){

	# Read DNAbin or file of gene.
	if(format == "phylip"){
		  if(aadata) data <- read.aa(sdata) else data <- read.dna(sdata)
	} else if(format == "fasta"){
	       	  if(aadata) data <- read.aa(sdata, format = "fasta") else data <- read.dna(sdata, format = "fasta")
	} else if(format == "bin"){
	          data <- sdata
	}


	# Run IQtree and extract the maximum likelihood, tree, and parameter estimates.
	
	if(!getTreeForced) ncore <- "AUTO"
	
	if(getTreeForced || !all(stats %in% c("chisq", "multlik", "biochemdiv", "maha"))) iqtreeres <- runIQtree(sdata, format = format, aadata = aadata, temp_name = geneName, iqtreePath = iqtreePath, model = model, tree = tree, ncore = ncore)
	
	# Get test statistics.
	
	results <- list()
	
	if("chisq" %in% stats){
         results$chisq <- get.chisqstat(data)
	}

    if("multlik" %in% stats){
         results$multlik <- get.unconstrained(data)
    }

    if("delta" %in% stats){
         results$delta <- iqtreeres$delta
    }

    if("biochemdiv" %in% stats){
         results$biocp <- get.biodivstat(data)
    }

    if("consind" %in% stats){
         if(!aadata) results$consind <- CI(iqtreeres$tree, phyDat(data)) else results$consind <- CI(iqtreeres$tree, phyDat(as.phyDat(data), type = "AA"))
    }

    if("brsup" %in% stats){
         results$brsup <- iqtreeres$meanNodeSupport
    }

    if("CIbrsup" %in% stats){
         results$CIbrsup <- iqtreeres$nodeSupport95
    }

    if("trlen" %in% stats){
         results$trlen <- iqtreeres$treeLength
    }
    
    if("imbal" %in% stats){
    	 results$imbal <- colless(as.treeshape(iqtreeres$tree))
    }

    if("stemmystat" %in% stats){
         results$stemmystat <- stemmy(iqtreeres$tree)
    }
    
    if("df" %in% stats){
    	 results$df <- get.df(iqtreeres$tree)
    }

    if("aindex" %in% stats){
    	 results$aindex <- iqtreeres$tree$edge.length / iqtreeres$treeLength
    }

    # Return test statistics, tree, and parameter estimates.
    
    if(getTreeForced || !all(stats %in% c("chisq", "multlik", "biochemdiv", "maha"))){
	
 	results$outputTree <- iqtreeres$tree
	
	if(length(grep("HKY|GTR", model)) == 1) results$piParams <- iqtreeres$piParams
	
	if(length(grep("[+]G", model)) == 1) results$alphaParam <- iqtreeres$alphaParam
	
	if(length(grep("GTR", model)) == 1) results$gtrMatrix <- iqtreeres$gtrMatrix

	if(length(grep("HKY", model)) == 1) results$trtvRatio <- iqtreeres$trtvRatio
	
    }
	
    return(results)

}