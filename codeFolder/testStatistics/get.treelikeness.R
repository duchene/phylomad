get.quartet.entropy <- function(pT1,pT2,pT3) -sum(c((pT1*log(pT1, 3)),(pT2*log(pT2, 3)), (pT3*log(pT3, 3))))

get.internode.cert <- function(pT1, pT2, pT3){
        dat <- sort(c(pT1, pT2, pT3), decreasing = T)
        res <- 1+((dat[1]/(dat[1]+dat[2]))*log2(dat[1]/(dat[1]+dat[2]))) + ((dat[2]/(dat[1]+dat[2]))*log2(dat[2]/(dat[1]+dat[2])))
        return(res)
}

## Add euclidean and chi-square!

get.dist2net <- function(pT1, pT2, pT3){
	concfactmat <- matrix(c(sort(c(pT1, pT2, pT3), decreasing = T), 0.499,0.499,0.002), 2, byrow=T)
	ait <- as.numeric(dist(concfactmat, method = "aitchison"))
	euc <- as.numeric(dist(concfactmat, method = "euclidean"))
	chi <- as.numeric(chisq.test(matrix(c(sort(c(pT1, pT2, pT3), decreasing = T), 0.5,0.5,0), 2, byrow=T))$statistic)
	dists2net <- c(aitchinson = ait, euclidean = euc, chisq = chi)
	return(dists2net)
}

get.dist2tr <- function(pT1, pT2, pT3){
	concfactmat <- matrix(c(sort(c(pT1, pT2, pT3), decreasing = T), 0.999,0.0005,0.0005), 2, byrow=T)
        ait <- as.numeric(dist(concfactmat, method = "aitchison"))
       	euc <- as.numeric(dist(concfactmat, method = "euclidean"))
        chi <- as.numeric(chisq.test(matrix(c(sort(c(pT1, pT2, pT3), decreasing = T), 1,0,0), 2, byrow=T))$statistic)
	dists2tr <- c(aitchinson = ait, euclidean = euc, chisq = chi)
        return(dists2tr)
}