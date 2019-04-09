print.bias.risk <- function(selectedStats, geneResults, Nsites, make.plot = F){

	outs <- vector()
	risktext <- vector()

	chisqthresholds <- data.frame(seqlen = c(500, 1000, 5000), minD = c(40, 70, 367), maxD = c(135, 258, 1273))
	multthresholds <- data.frame(seqlen = c(200, 1000, 5000), minD = c(42.57623, 80.66334, 141.216), maxD = c(116.3087, 190.2008, 247.9695))
	biochthresholds <- data.frame(seqlen = c(200, 1000, 5000), minD = c(-0.9583055, -2.3559903, -6.14660146), maxD = c(-4.262108, -8.001475, -14.573290))
	consthresholds <- data.frame(seqlen = c(200, 1000, 5000), minD = c(-0.6097077, -1.250933, -2.509672), maxD = c(-2.21697, -3.555679, -5.434672))
	mahathresholds <- data.frame(seqlen = c(200, 1000, 5000), minD = c(6.012026, 17.662311, 26.279046))
	
	basicthresholds <- do.call("cbind", list(chisqthresholds, multthresholds, biochthresholds[,2:3], consthresholds[,2:3], mahathresholds[,2]))
	colnames(basicthresholds) <- c("seqlen.chisq", "minD.chisq", "maxD.chisq", "seqlen", "minD.multlik", "maxD.multlik", "minD.biochemdiv", "maxD.biochemdiv", "minD.consind", "maxD.consind", "minD.maha")
	
	thresholds <- matrix(NA, nrow = 2, ncol = length(selectedStats))
	rownames(thresholds) <- c("mid.risk.threshold", "high.risk.threshold")
	colnames(thresholds) <- selectedStats

	if("chisq" %in% selectedStats){
               thresholdMidRisk <- predict(lm(minD ~ seqlen, data = chisqthresholds), data.frame(seqlen = Nsites))
               thresholdHighRisk <- predict(lm(maxD ~ seqlen, data = chisqthresholds), data.frame(seqlen = Nsites))
               if(is.na(geneResults$chisq.sdpd)){
			risktext <- c(risktext, "This locus is at UNKNOWN risk of leading to biased inferences due to compositional heterogeneity. It looks like the assessment failed.", "")
                        outs <- c(outs, "Unknown.risk")
	       } else if(geneResults$chisq.sdpd > thresholdHighRisk){
                        risktext <- c(risktext, "This locus is at HIGH risk of leading to biased inferences due to compositional heterogeneity. It is highly unadvisable to analyse this locus using homogeneous substitution models such as those of the GTR family.", "")
                        outs <- c(outs, "high.risk")
               } else if(geneResults$chisq.sdpd > thresholdMidRisk){
                        risktext <- c(risktext, "This locus is at MEDIUM risk of leading to biased inferences due to compositional heterogeneity. These data might not lead to biased inferences when using homogeneous substitution models such as those of the GTR family. However, these inferences should be made with caution and the results of this assessment should be reported. Alternatively, a heterogeneous substitution model can be used.", "")
                        outs <- c(outs, "mid.risk")
               } else {
                        risktext <- c(risktext, "This locus is at LOW risk of leading to biased inferences due to compositional heterogeneity. Homogeneous substitution models such as those of the GTR family might not lead to biased inferences due to compositional heterogeneity.", "")
                        outs <- c(outs, "low.risk")
               }
               names(outs)[length(outs)] <- "chi.squared.assessment"
	       thresholds[1,"chisq"] <- thresholdMidRisk
	       thresholds[2,"chisq"] <-	thresholdHighRisk	       
    }

	if("multlik" %in% selectedStats){
	       thresholdMidRisk <- predict(lm(minD ~ log(seqlen), data = multthresholds), list(seqlen = Nsites))
	       thresholdHighRisk <- predict(lm(maxD ~ log(seqlen), data = multthresholds), list(seqlen = Nsites))
	       if(is.na(geneResults$multlik.sdpd)){
			risktext <- c(risktext, "This locus is at UNKNOWN risk of leading to biased inferences according to assessment using the multinomial likelihood. It looks like the assessment failed.", "")
                        outs <- c(outs, "Unknown.risk")
	       } else if(geneResults$multlik.sdpd > thresholdHighRisk){
	       		risktext <- c(risktext, "This locus is at HIGH risk of leading to based inferences according to assessment using the multinomial likelihood. The model is underparameterized.", "")
			outs <- c(outs, "high.risk")
	       } else if(geneResults$multlik.sdpd > thresholdMidRisk){
			risktext <- c(risktext, "This locus is at MEDIUM risk of leading to based inferences according to assessment using the multinomial likelihood. This indicates that the model is underparameterized, but might still provide some reasonable inferences.", "")
                        outs <- c(outs, "mid.risk")
	       } else {
			risktext <- c(risktext, "This locus is at LOW risk of leading to based inferences according to assessment using the multinomial likelihood. This suggests that the model is not underparameterized.", "")
                        outs <- c(outs, "low.risk")
	       }
	       names(outs)[length(outs)] <- "multinomial.likelihood.assessment"
	       		thresholds[1,"multlik"] <- thresholdMidRisk
		        thresholds[2,"multlik"] <- thresholdHighRisk
	}

	if("biochemdiv" %in% selectedStats){
	       thresholdMidRisk <- predict(lm(minD ~ log(seqlen), data = biochthresholds), list(seqlen = Nsites))
               thresholdHighRisk <- predict(lm(maxD ~ log(seqlen), data = biochthresholds), list(seqlen = Nsites))
	       if(is.na(geneResults$biochemdiv.sdpd)){
			risktext <- c(risktext, "This locus is at UNKNOWN risk of leading to biased inferences according to the biochemical diversity. It looks like the assessment failed.", "")
                        outs <- c(outs, "Unknown.risk")
	       } else if(geneResults$biochemdiv.sdpd < thresholdHighRisk){
                        risktext <- c(risktext, "This locus is at HIGH risk of leading to biased inferences according to the biochemical diversity. Some explanations include excessive substitutional saturation or a highly variable subsitution process across branches and across sites. This test might also reject the model if the locus is extremely short.", "")
                        outs <- c(outs, "high.risk")
               } else if(geneResults$biochemdiv.sdpd < thresholdMidRisk){
                        risktext <- c(risktext, "This locus is at MEDIUM risk of leading to biased inferences according to the biochemical diversity. Some explanations include considerable substitutional saturation or a variable subsitution process across branches and across sites.", "")
                        outs <- c(outs, "mid.risk")
               } else {
                        risktext <- c(risktext, "This locus is at LOW risk of leading to biased inferences according to the biochemical diversity. This suggests that there might not be excessive saturation or the substitution process might be relatively constant across branches and across sites.", "")
                        outs <- c(outs, "low.risk")
               }
	       names(outs)[length(outs)] <- "biochemical.diversity.assessment"
	       thresholds[1,"biochemdiv"] <- thresholdMidRisk
	       thresholds[2,"biochemdiv"] <- thresholdHighRisk
	}

	if("consind" %in% selectedStats){
	       thresholdMidRisk <- predict(lm(minD ~ log(seqlen), data = consthresholds), list(seqlen = Nsites))
	       thresholdHighRisk <- predict(lm(maxD ~ log(seqlen), data = consthresholds), list(seqlen = Nsites))
	       if(is.na(geneResults$consind.sdpd)){
			risktext <- c(risktext, "This locus is at UNKNOWN risk of leading to biased inferences according to the consistency index. It looks like the assessment failed.", "")
                        outs <- c(outs, "Unknown.risk")
               } else if(geneResults$consind.sdpd < thresholdHighRisk){
                        risktext <- c(risktext, "This locus is at HIGH risk of leading to biased inferences according to the consistency index. One explanation is that there is strong variation in the rates of evolution across both branches and sites. This test might also reject the model if the locus is extremely short.", "")
                        outs <- c(outs, "high.risk")
               } else if(geneResults$consind.sdpd < thresholdMidRisk){
                        risktext <- c(risktext, "This locus is at MEDIUM risk of leading to biased inferences according to the consistency index. One explanation is that there is considerable variation in the rates of evolution across both branches and sites.", "")
                        outs <- c(outs, "mid.risk")
               } else {
                        risktext <- c(risktext, "This locus is at LOW risk of leading to biased inferences according to the consistency index. This suggests that the substitution process might be relatively constant across branches and across sites.", "")
                        outs <- c(outs, "low.risk")
               }
	       names(outs)[length(outs)] <- "consistency.index.assessment"
	       thresholds[1,"consind"] <- thresholdMidRisk
               thresholds[2,"consind"] <- thresholdHighRisk
	}

	if("maha" %in% selectedStats){
	       thresholdMidRisk <- predict(lm(minD ~ log(seqlen), data = mahathresholds), list(seqlen = Nsites))
	       if(is.na(geneResults$maha.sdpd)){
			risktext <- c(risktext, "This locus is at UNKNOWN risk of leading to biased inferences according to the mahalanobis summary statistic. It looks like some part of the assessment failed.", "")
                        outs <- c(outs, "Unknown.risk")
	       } else if(geneResults$maha.sdpd > thresholdMidRisk){
                        risktext <- c(risktext, "This locus is at MODERATE risk of leading to biased inferences according to the mahalanobis distance. This means that one of the other test statistics considered suggests a substantial amount of bias, or that multiple test statistics suggest at least some risk of biased inferences.", "")
                        outs <- c(outs, "mid.risk")
               } else {
                        risktext <- c(risktext, "This locus is at LOW risk of leading to biased inferences according to the mahalanobis summary statistic. This suggests that the model unlikely to lead to biased inferences according to the group of test statistics considered.", "")
                        outs <- c(outs, "low.risk")
               }
	       names(outs)[length(outs)] <- "mahalanobis.assessment"
	       thresholds[1,"maha"] <- thresholdMidRisk
	}
	
	risktext <- c(risktext, "This advice is based on simulations in the following studies:", "Duchêne, D.A., Duchêne, S., & Ho, S.Y.W. (2017). New Statistical Criteria Detect Phylogenetic Bias Caused by Compositional Heterogeneity. Molecular Biology and Evolution, 34(6), 1529-1534.", "Duchêne, D.A., Duchêne, S., & Ho, S.Y.W. (in prep). New Statistical Criteria Detect Biased Inferences From Phylogenomic Data. ")
	
	thresholds <- thresholds[, selectedStats]
	writeLines(risktext, con = "bias.risk.report.txt")
	outs <- list(reccommendation = matrix(c(thresholds[1,], thresholds[2,], outs), nrow = 1, dimnames = list("", c(paste0(selectedStats, ".mid.risk.threshold"), paste0(selectedStats, ".high.risk.threshold"), names(outs)))), thresholds = basicthresholds)
	return(outs)

}
