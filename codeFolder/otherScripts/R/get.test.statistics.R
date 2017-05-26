get.test.statistics <- function(sdata, format = "phylip", aadata = F, geneName = "empirical", phymlPath, model = "GTR+G", stats = c("chisq", "multlik", "delta", "biochemdiv", "consind", "brsup", "trlen", "maha"), tree = NULL){

	# Read DNAbin or file of gene.
	if(format == "phylip"){
		  if(aadata) data <- read.aa(sdata) else data <- read.dna(sdata)
	} else if(format == "fasta"){
	       	  if(aadata) data <- read.aa(sdata, format = "fasta") else data <- read.dna(sdata, format = "fasta")
	} else if(format == "bin"){
	          data <- sdata
	}


	# Run PhyML and extract the maximum likelihood, tree, and parameter estimates.
	
	phymlres <- runPhyML(sdata, format = format, aadata = aadata, temp_name = geneName, phymlPath = phymlPath, model = model, tree = tree)
	
	# Get test statistics.
	
	results <- list()
	
	if("chisq" %in% stats){
         results$chisq <- get.chisqstat(data)
	}

    if("multlik" %in% stats){
         results$multlik <- phymlres$uncLikelihood
    }

    if("delta" %in% stats){
         results$delta <- phymlres$delta
    }

    if("biochemdiv" %in% stats){
         results$biocp <- get.biodivstat(data)
    }

    if("consind" %in% stats){
         if(!aadata) results$consind <- CI(phymlres$tree, phyDat(data)) else results$consind <- CI(phymlres$tree, phyDat(as.phyDat(data), type = "AA"))
    }

    if("brsup" %in% stats){
         results$brsup <- phymlres$meanNodeSupport
    }

    if("CIbrsup" %in% stats){
         results$CIbrsup <- phymlres$nodeSupport95
    }

    if("trlen" %in% stats){
         results$trlen <- phymlres$treeLength
    }
    
    if("imbal" %in% stats){
    	 results$imbal <- colless(as.treeshape(phymlres$tree))
    }

    if("stemmystat" %in% stats){
         results$stemmystat <- stemmy(phymlres$tree)
    }
    
    if("df" %in% stats){
    	 results$df <- get.df(phymlres$tree)
    }

    if("aindex" %in% stats){
    	 results$aindex <- phymlres$tree$edge.length / phymlres$treeLength
    }

    # Return test statistics, tree, and parameter estimates.
	
 	results$outputTree <- phymlres$tree
	
	if(length(grep("HKY|GTR", model)) == 1) results$piParams <- phymlres$piParams
	
	if(length(grep("[+]G", model)) == 1) results$alphaParam <- phymlres$alphaParam
	
	if(length(grep("GTR", model)) == 1) results$gtrMatrix <- phymlres$gtrMatrix

	if(length(grep("HKY", model)) == 1) results$trtvRatio <- phymlres$trtvRatio
	
	return(results)


}