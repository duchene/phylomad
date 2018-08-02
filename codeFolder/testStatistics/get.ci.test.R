require(phangorn)

get.ci.test <- function(phy, al){
	
	invarsites <- which(apply(as.character(al), 2, function(x) length(unique(x))) == 1)
	if(length(invarsites) > 0) al <- al[,-invarsites] 

	ci <- CI(phy, phyDat(al), sitewise = T)^-1

	m <- phangorn:::lowerBound(subset(phyDat(al), phy$tip.label))
	weight <- attr(phyDat(al), "weight")
	s <- sum(base.freq(al)**2)
        k <- nrow(al)
        null_changes <- (1-s)*k
	Eci <- null_changes / mean(m)

	test <- t.test(ci, mu = Eci)

	res <- list()
	res$t <- test$statistic
	res$p <- test$p.value
	return(res)
}