runPhyML <- function(sdata, format = 'phyllip', temp_name, phymlPath = '~/Downloads/PhyML-3.1/PhyML-3.1_macOS-MountainLion', model = 'JC'){
    require(phangorn)

    if(format == 'fasta'){
        d <- read.dna(sdata, format = 'fasta')
        fileName = gsub('fasta', 'phy', sdata)
        write.dna(d, file = fileName)
    }else if(format == 'DNAbin'){
        write.dna(sdata, temp_name)
        fileName = temp_name
    }else{
        fileName = sdata
    }
    print(fileName)
    if(model == 'JC'){
        phymlOptions = ' -m jc69 -c 1 --q -i '
    }
    else if(model == 'GTR+G'){
        phymlOptions = ' -m gtr -a e --q -i '
    }

    phymlCommand = paste0(phymlPath, phymlOptions, fileName)

    system(phymlCommand)
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

    if(model == 'GTR+G'){
        gtrLocation <- grep('GTR relative rate parameters', outputStats)
        gtrMatrix <- as.numeric(gsub('[A-Z]| |[-]|>|<', '',  outputStats[(gtrLocation+1):(gtrLocation+6)]))

        piLocation <- grep('Nucleotides frequencies', outputStats)
        piParams <- as.numeric(gsub('[A-Z]|[a-z]| |[-]|>|<|[(]|[)]|=', '', outputStats[(piLocation+1):(piLocation+4)]))

        alphaLocation <- grep('Gamma shape parameter', outputStats)
        alphaParam <- as.numeric(gsub('[A-Z]|[a-z]| |:|[-]|\t', '', outputStats[alphaLocation]))

        testStats$gtrMatrix <-  gtrMatrix
        testStats$piParams <-  piParams
        testStats$alphaParam <- alphaParam
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
