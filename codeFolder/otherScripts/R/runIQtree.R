# This function extracts some of the main results from an IQtree run output file.

runIQtree <- function(sdata, format = 'phylip', aadata = aadata, temp_name, iqtreePath = '~/Downloads/iqtree-1.7-beta9-MacOSX/bin/iqtree', model = 'JC', tree = NULL, ncore = 1, symtest = F){
       
       require(phangorn)

    if(format == 'fasta'){
        if(aadata) d <- read.aa(sdata, format = 'fasta') else d <- read.dna(sdata, format = 'fasta')
        fileName <- gsub('fasta', 'phy', sdata)
        write.dna(d, file = fileName)
    } else if(format == 'bin'){
        write.dna(sdata, file = temp_name)
        fileName <- temp_name
    } else if(format == "nexus"){
        if(aadata) d <- as.AAbin(read.nexus.data(sdata)) else d <- as.DNAbin(read.nexus.data(sdata))
        fileName <- gsub('nexus', 'phy', sdata)
        write.dna(d, file = fileName)
    } else {
        fileName <- sdata
    }
    print(paste("Locus", fileName, "was read successfully."))

    if(is.null(tree)){
        iqtreeCommand = paste0(iqtreePath, if(!aadata) " -st DNA " else " -st AA ", " -s ", fileName, " -m ", model, " -alrt 1000 -fast -nt ", ncore, if(symtest) " --symtest")
    } else {
        treenumber <- round(runif(1, min = 1000, max = 9999))
	treefile <- paste0("temp.", treenumber, ".tre")
        write.tree(tree, file = treefile)
        iqtreeCommand = paste0(iqtreePath, if(!aadata) " -st DNA " else " -st AA ", " -s ", fileName, " -m ", model, " -g ", treefile, " -alrt 1000 -nt ", ncore, if(symtest) " --symtest")
    }

    system(iqtreeCommand)
    if(!is.null(tree)) system(paste0("rm ", treefile))
    
    allout <- readLines(paste0(fileName, ".iqtree"))

	torm <- grep(paste0(fileName, c(".ckp.gz", ".bionj", ".log", ".mldist", ".treefile", ".uniqueseq.", ".contree", ".splits.nex", ".iqtree", ".parstree"), collapse = "|"), dir(), value = T)
	
	for(i in torm) system(paste0("rm ", i))
       
       res <- list()
       
       esttree <- read.tree(text = allout[grep("Tree in newick format:", allout)[1] + 2])

       boots <- as.numeric(esttree$node.label)
       res$meanNodeSupport <- mean(boots, na.rm = T)
       
       res$nodeSupport95 <- diff(quantile(boots, c(0.025, 0.975), na.rm = T))

       trlen <- grep("Total tree length", allout, value = T)
       res$treeLength <- as.numeric(gsub(".* ", "", trlen))

       if(length(grep("HKY|GTR", model))==1){
        	piLocation <- grep('State frequencies', allout)+2
        	piParams <- as.numeric(gsub('.* ', '', allout[(piLocation):(piLocation+3)]))
        	res$piParams <-  piParams
       }

       if(length(grep("[+]G", model))==1){
        	alphaLocation <- grep('Gamma shape alpha', allout)
        	alphaParam <- as.numeric(gsub(".* ", "", allout[alphaLocation]))
       		res$alphaParam <- alphaParam
       }

       if(length(grep("GTR", model))==1){
        	gtrLocation <- grep('Rate parameter R', allout)+2
        	gtrMatrix <- as.numeric(gsub('.* ', '',  allout[(gtrLocation):(gtrLocation+5)]))
        	res$gtrMatrix <-  gtrMatrix
       }

       if(length(grep("HKY", model))==1){
        	trtvLocation <- grep('A-G:', allout)
        	trtvRatio <- as.numeric(gsub(".* ", "", allout[trtvLocation]))
        	res$trtvRatio <-  trtvRatio
       }

       likse <- grep("Log-likelihood of the tree:", allout, value = T)
       res$likelihood <- as.numeric(gsub(".*:|[(].*| ", "", likse))
       
       unclik <- grep("Unconstrained log-likelihood", allout, value = T) 
       res$uncLikelihood <- as.numeric(gsub(".* ", "", unclik))
       
       res$delta <- res$uncLikelihood - res$likelihood
       
       res$tree <- esttree
       
       # Other parameters that might be of interest in the future
       #params <- grep("Number of free parameters", allout, value  = T)
       #res$Nparams <- gsub(".* ", "", params)       
       #aic <- grep("Akaike information criterion", allout, value = T)
       #res$aic <- gsub(".* ", "", aic[2])
       #res$aicc <- gsub(".* ", "", aic[3])       
       #bic <- grep("Bayesian information criterion", allout, value = T)
       #res$bic <- gsub(".* ", "", bic[2])       
       #bootrange <- range(boots, na.rm = T)
       #res$bootmax <- max(bootrange)
       #res$bootmin <- min(bootrange)
       #res$bootrangesize <- diff(bootrange)       
       #res$likelihoodSE <- gsub(".*[(]|[)].*|.* ", "", likse)

       return(res)
}