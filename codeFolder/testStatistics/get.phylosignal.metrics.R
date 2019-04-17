# The input of these functions is a vector of quartet probabilities, qp, e.g., c(pT1, pT2, pT3). The last function, get.binom.p takes the same vector but with the number of trials appended as the last value.

get.quartet.entropy <- function(qp) -sum(c((qp[1]*log(qp[1], 3)),(qp[2]*log(qp[2], 3)), (qp[3]*log(qp[3], 3))))

get.internode.cert <- function(qp){
        dat <- sort(qp, decreasing = T)
        res <- 1+((dat[1]/(dat[1]+dat[2]))*log2(dat[1]/(dat[1]+dat[2]))) + ((dat[2]/(dat[1]+dat[2]))*log2(dat[2]/(dat[1]+dat[2])))
        return(res)
}

get.dist2net <- function(qp){
	concfactmat <- matrix(c(sort(qp, decreasing = T), 0.499,0.499,0.002), 2, byrow=T)
	#ait <- as.numeric(dist(concfactmat, method = "aitchison"))
	euc <- as.numeric(dist(concfactmat, method = "euclidean"))
	chi <- as.numeric(suppressWarnings(chisq.test(matrix(c(sort(qp, decreasing = T), 0.5,0.5,0), 2, byrow=T)))$statistic)
	dists2net <- c(aitchinson = NA, euclidean = euc, chisq = chi)
	return(dists2net)
}

get.dist2tr <- function(qp){
	concfactmat <- matrix(c(sort(qp, decreasing = T), 0.999,0.0005,0.0005), 2, byrow=T)
        #ait <- as.numeric(dist(concfactmat, method = "aitchison"))
       	euc <- as.numeric(dist(concfactmat, method = "euclidean"))
        chi <- as.numeric(suppressWarnings(chisq.test(matrix(c(sort(qp, decreasing = T), 1,0,0), 2, byrow=T)))$statistic)
	dists2tr <- c(aitchinson = NA, euclidean = euc, chisq = chi)
	return(dists2tr)
}

get.binom.p <- function(qpNtrials){
	bipartsuccs <- sort(qpNtrials[1:3]*qpNtrials[4], decreasing = T)
	pval <- binom.test(round(bipartsuccs[1]), round(qpNtrials[4]), 1/3, alternative = "greater")$p.value
	return(pval)
}

get.dstat <- function(qp){
        dat <- sort(qp, decreasing = T)
        if(dat[2] == dat[3]) return(NA)
        dstat <- (dat[2] - dat[3]) / (dat[2] + dat[3])
        return(dstat)
}

get.kcstat <- function(qp){
        dat <- sort(qp, decreasing = T)
        if(dat[2] == dat[3]) return(NA)
        kcstat <- (dat[1] - dat[3]) / (dat[2] - dat[3])
        return(kcstat)
}