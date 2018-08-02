require(phangorn)

get.comp.test <- function(phy, al){
	
	invarsites <- which(apply(as.character(al), 2, function(x) length(unique(x))) == 1)
        if(length(invarsites) > 0) al <- al[,-invarsites]
	
	phy$node.label <- NULL
	phy$edge.length <- NULL
	taxord <- gsub("[(]|[)]|;", "", write.tree(phy))
	taxord <- strsplit(taxord, ",")[[1]]
	order_alignment <- al[taxord,]
	mat <- as.matrix(al)
	sitewise_rle = apply(mat, 2, function(x) length(rle(x)$lengths))

	total_rle <- sum(sitewise_rle)
	mean_rle <- mean(sitewise_rle)
	
	s <- sum(base.freq(al)**2)
	k <- nrow(al)
	p <- 1 - pbinom(mean_rle, k, 1-s)

	null_mean_rle <- (1-s)*k + 1
	test <- t.test(sitewise_rle, mu=null_mean_rle)
	
	res <- list()
	res$t <- test$statistic
	res$p <- test$p.value
	
	return(res) 
}