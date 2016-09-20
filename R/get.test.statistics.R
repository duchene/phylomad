get.test.statistics <- function(sdata, format = "phyllip", geneName = "empirical", phymlPath, model = "GTR+G"){

	# Read DNAbin or file of gene.
	if(format == "phyllip"){
		  data <- read.dna(sdata)
	} else if(format == "fasta"){
	       	  data <- read.dna(sdata, format = "fasta")
	} else if(format == "DNAbin"){
	          data <- sdata
	}


	# Run PhyML and extract the maximum likelihood, tree, and parameter estimates.
	
	phymlres <- runPhyML(sdata, format = format, temp_name = geneName, phymlPath = phymlPath, model = model)
	
	# Get test statistics.
	chisq <- getchisqs(data)
	
	biochem <- biochemdiv.stat(data)

	# Return test statistics, tree, and parameter estimates.
	
	results <- list(multinomial.lik = phymlres$uncLikelihood, chisq.stat = chisq, biochemdiv = biochem, mean.branch.sup = phymlres$meanNodeSupport, CI.branch.sup = phymlres$nodeSupport95, delta = phymlres$delta, tree.length = phymlres$treeLength, outputTree = phymlres$tree)
	
	if(model == "GTR+G"){
		 results$gtrMatrix <- phymlres$gtrMatrix
		 results$piParams <- phymlres$piParams
		 results$alphaParam <- phymlres$alphaParam 
	}
	
	return(results)


}