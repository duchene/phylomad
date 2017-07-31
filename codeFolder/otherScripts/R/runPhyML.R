runPhyML <- function(sdata, format = 'phylip', aadata = aadata, temp_name, phymlPath = '~/Downloads/PhyML-3.1/PhyML-3.1_macOS-MountainLion', model = 'JC', tree = NULL){
    require(phangorn)
    
    if(format == 'fasta'){
    	if(aadata) d <- read.aa(sdata, format = 'fasta') else d <- read.dna(sdata, format = 'fasta')
        fileName = gsub('fasta', 'phy', sdata)
	write.dna(d, file = fileName)
    } else if(format == 'bin'){
        write.dna(sdata, temp_name)
        fileName = temp_name
    } else if(format == "nexus"){
        if(aadata) d <- as.AAbin(read.nexus.data(sdata)) else d <- as.DNAbin(read.nexus.data(sdata))
	fileName = gsub('nexus', 'phy', sdata)
	write.dna(d, file = fileName)
    } else {
        fileName = sdata
    }
    print(paste("Locus", fileName, "was read successfully."))

    if(length(grep("[+]G", model)) == 1) RAS <- " -a e " else RAS <- " -c 1 "

    if(length(grep('JC', model)) == 1){
        phymlOptions = 'jc69'
    } else if(length(grep('HKY', model)) == 1){
       phymlOptions = 'hky85'
    } else if(length(grep('GTR', model)) == 1){
        phymlOptions = 'gtr'
    } else if(length(grep('LG', model)) == 1){
        phymlOptions = 'lg'
    } else if(length(grep('WAG', model)) == 1){
        phymlOptions = 'wag'
    } else if(length(grep('JTT', model)) == 1){
        phymlOptions = 'jtt'
    } else if(length(grep('Dayhoff', model)) == 1){
        phymlOptions = 'dayhoff'
    }

    if(aadata) aasetting <- " -d aa" else aasetting <- ""

    if(is.null(tree)){
    	phymlCommand = paste0(phymlPath, aasetting, " -m ", phymlOptions, RAS, "--q -i ", fileName)
    } else {
      	treenumber <- round(runif(1, min = 1000, max = 9999))
	write.tree(tree, file = paste0("temp.", treenumber, ".tre"))
        phymlCommand = paste0(phymlPath, aasetting, " -m ", phymlOptions, " -o lr -u temp.", treenumber, ".tre", RAS, "--q -i ", fileName)
    }

    system(phymlCommand)
    if(!is.null(tree)) system(paste0("rm temp.", treenumber, ".tre"))
    outTreeName <- paste0(fileName, '_phyml_tree.txt')
    outStatsName <- paste0(fileName, '_phyml_stats.txt')
    outputStats <- readLines(outStatsName)
    outputTree <- read.tree(outTreeName)

    testStats <- list()

    meanNodeSupport <- mean(as.numeric(outputTree$node.label), na.rm = T)
    nodeSupport95 <- diff(quantile(as.numeric(outputTree$node.label), c(0.025, 0.975), na.rm = T))
    treeLength <- sum(outputTree$edge.length)

    testStats$meanNodeSupport <-  meanNodeSupport
    testStats$nodeSupport95 <-  nodeSupport95
    testStats$treeLength <- treeLength

    if(length(grep("HKY|GTR", model))==1){
    	piLocation <- grep('Nucleotides frequencies', outputStats)
        piParams <- as.numeric(gsub('[A-Z]|[a-z]| |[-]|>|<|[(]|[)]|=', '', outputStats[(piLocation+1):(piLocation+4)]))
	testStats$piParams <-  piParams
    }
    
    if(length(grep("[+]G", model))==1){
    	alphaLocation <- grep('Gamma shape parameter', outputStats)
        alphaParam <- as.numeric(gsub('[A-Z]|[a-z]| |:|[-]|\t', '', outputStats[alphaLocation]))
	testStats$alphaParam <- alphaParam
    }
    
    if(length(grep("GTR", model))==1){
    	gtrLocation <- grep('GTR relative rate parameters', outputStats)
        gtrMatrix <- as.numeric(gsub('[A-Z]| |[-]|>|<', '',  outputStats[(gtrLocation+1):(gtrLocation+6)]))
	testStats$gtrMatrix <-  gtrMatrix
    }

    if(length(grep("HKY", model))==1){
        trtvLocation <- grep('Transition', outputStats)
	trtvRatio <- as.numeric(gsub('[A-Z]|[a-z]| |:|[-]|\t|[/]|^[.]', '',  outputStats[trtvLocation]))
        testStats$trtvRatio <-  trtvRatio
    }    

    likelihoodLocation <- grep('Log-likelihood', outputStats)
    likelihood <- as.numeric(gsub('[A-Z]|[a-z]| |:|\t|[/]|[.] |[a-z][-][a-z]', '', outputStats[likelihoodLocation]))

    uncLikelihoodLocation <- grep('Unconstrained likelihood', outputStats)
    uncLikelihood <- as.numeric(gsub('[A-Z]|[a-z]| |:|[.] |\t', '', outputStats[uncLikelihoodLocation]))

    delta <- uncLikelihood - likelihood 

    testStats$likelihood <-  likelihood
    testStats$uncLikelihood <-  uncLikelihood
    testStats$delta <-  delta
    testStats$tree <- outputTree

    system(paste('rm', outTreeName))
    system(paste('rm', outStatsName))

    return(testStats)
}
