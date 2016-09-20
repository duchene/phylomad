# This function runs model adequacy assessment for genes in a file.


run.genome <- function(format = "phyllip", model = "GTR+G", phymlPath, Nsims = 100, paragenome = F, paragene = F, ncore = 1){
	   require(phangorn)
	   
	   if(format == "phyllip"){
	   	   genes <- grep(".+phy", dir(), value = T)
		   
	   } else if(format == "fasta"){
	     	   genes <- grep(".+[.]fasta", dir(), value = T)
	   }
	   geneStats <- list()
	   
	   if(!paragenome){
	     for(i in 1:length(genes)){
	   	 print(genes[i])
		 genes[i]
	   	 geneStats[[i]] <- run.gene(genes[i], format = format, model = model, phymlPath = phymlPath, Nsims = Nsims, para = paragene, ncore = ncore)
	     }
	   }else{
	     ## START PARALLEL COMPUTING
	     runGene <- function(i){
	       return(run.gene(genes[i], format = format, model = model, phymlPath = phymlPath, Nsims = Nsims))
	     }
	     
	     cl <- makeCluster(ncore)
	     registerDoParallel(cl)
	     geneRuns <- foreach(x = 1:length(genes), .packages = c('phangorn', 'ape'), .export = c('get.test.statistics', 'runPhyML', 'getchisqs', 'run.gene')) %dopar% runGene(x)
	     stopCluster(cl)
	     geneStats <- geneRuns
	     ## END PARALLEL COMPUTING
	   }
	   
	   names(geneStats) <- genes

	   genomeStats <- matrix(NA, nrow = length(genes), ncol = 14)
	   for(i in 1:14){
	   	 genomeStats[,i] <- sapply(geneStats, function(x) x[[i]])
 	   }
	   geneTrees <- lapply(geneStats, function(x) x[[15]])
	   names(geneTrees) <- genes
	   
	   if(model == "GTR+G"){
	   	    geneInferences <- list()
		    for(i in 1:length(geneStats)){
		    	  print(geneStats[[i]][16:18])
		    	  geneInferences[[i]] <- geneStats[[i]][16:18]
		    }
		    names(geneInferences) <- genes
	   }	   


	   rownames(genomeStats) <- genes
	   colnames(genomeStats) <- names(geneStats[[1]])[1:14]
	   results <- list(genome.results = genomeStats, empirical.trees = geneTrees, empirical.parameters = geneInferences)

	   return(results)

}