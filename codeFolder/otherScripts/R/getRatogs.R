# This function is to get ratograms from the output of beast2. It takes an input file, and returns an R object of trees.

require(phangorn)

getRatogs <- function(treesf){
	  trs <- readLines(treesf)
	  trs <- gsub("[[]&([a-z])+[=]", ":", trs)
	  trs <- gsub("[]]:([0-9]|[.])+", "", trs)
	  
	  #trs <- writeLines(trs)
	  #trs <- read.nexus(trs)

	  trsname <- paste0(sample(letters, 5), collapse = "")
	  #print(trsname)
	  writeLines(trs, con = trsname)
	  trs <- read.nexus(trsname)
	  system(paste("rm", trsname))
	  return(trs)
}