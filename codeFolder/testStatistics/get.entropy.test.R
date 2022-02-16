require(ape)


get.entropy.test <- function(al, only.varsites = F){

### Create functions needed ###

calc_term <- function(n, p_i){

  res = c()
  for(i in seq(1,n,1)){
    res[i] = choose(n,i)*p_i^(i)*((1-p_i)^(n -i))*lfactorial(i)
  }
  
  return(sum(res))
}

multinomial_entropy <- function(n, p){
  
  term_3 = sum(sapply(p, FUN = calc_term, n = n))
  calc <- -lfactorial(n) - n*sum(p*log(p)) + term_3
    
  return(calc)
}

smart_table <- function(x){
  bases <- c("a", "c", "g", "t")
  return(table(factor(x, levels = bases)))
}

multinomial_info <- function(x, p){ 
  P_x = dmultinom(x, prob = p)
  return(-log(P_x))
}

calculate_sitewise_info <- function(mat, base_freqs){
  site_wise_counts <- apply(mat, 2, FUN = smart_table)
  site_wise_info <- apply(site_wise_counts, 2, multinomial_info, p = base_freqs)
  return(site_wise_info)
}

calculate_estimated_entropy <- function(site_wise_info, mat){
  estimated_entropy <- sum(site_wise_info)/dim(mat)[2]

  return(estimated_entropy)
}

### Functions needed created. Statstic calculated below ###
  
  #if(only.varsites) al <- al[,seg.sites(al)]
  if(only.varsites) al <- al[,pis(al, "index")]
  if(ncol(al) < 20){
        res <- list(estent = NA, nullent = NA, p = NA, tstat = NA)
        return(res)
  }
  mat <- tolower(as.character(al))
  
  base_freqs <- base.freq(al) #observed base frequencies #
  
  if(any(base_freqs > 0.5)){
        res <- list(estent = NA, nullent = NA, p = NA, tstat = NA)
        return(res)
  }
  
  n <- dim(al)[1] # the number of taxa# 
  num_seqs <- dim(al)[2] #number of sites#
  
  null_entropy <- multinomial_entropy(n, base_freqs)  #the entropy expected under saturation#
  
  site_wise_info <- calculate_sitewise_info(mat, base_freqs)
  
  estimated_entropy <- calculate_estimated_entropy(site_wise_info, mat)
  
  if(var(site_wise_info) == 0){
        if(site_wise_info[1] > null_entropy) tt <- list(statistic = NA, p.value = 0)
  } else {
    	tt <- t.test(site_wise_info, mu = null_entropy, alternative = "two.sided")
  }
  
  estent <- estimated_entropy
  nullent <- null_entropy
  
  res <- list()
  res$t <- tt$statistic
  res$p <- tt$p.value
  return(res)
}
