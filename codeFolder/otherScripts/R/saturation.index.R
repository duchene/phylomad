require(ape)

calc_term <- function(n, p_i){

  res = c()
  for(i in seq(1,n,1)){
    res[i] = choose(n,i)*p_i^(i)*((1-p_i)^(n -i))*log(factorial(i))
  }
  
  return(sum(res))
}

multinomial_entropy <- function(n, p){
  
  term_3 = sum(sapply(p, FUN = calc_term, n = n))
  calc <- -log(factorial(n)) - n*sum(p*log(p)) + term_3
    
  return(calc)
}

smart_table <- function(x){
  bases <- c("a", "c", "g", "t")
  return(table(factor(x, levels = bases)))
}


calculate_sitewise_info <- function(mat, base_freqs){
  site_wise_counts <- apply(mat, 2, FUN = smart_table)
  print(site_wise_counts[0])

  site_wise_info <- apply(site_wise_counts, 2, multinomial_info, p = base_freqs)
  
  return(site_wise_info)
}

calculate_estimated_entropy <- function(site_wise_info){
  estimated_entropy <- sum(site_wise_info)/dim(mat)[2]
  
  return(estimated_entropy)
  
}

calculate_index <- function(al, give_p = TRUE){
  mat <- as.character(al)
  
  base_counts <- as.vector(table(mat))
  base_freqs <- base_counts / prod(dim(mat))
  base_freqs #observed base frequencies #
  
  n <- dim(mat)[1] # the number of taxa# 
  num_seqs <- dim(mat)[2] #number of sites#
  
  null_entropy <- multinomial_entropy(n, base_freqs)  #the entropy expected under saturation#
  
  site_wise_info <- calculate_sitewise_info(mat, base_freqs)
  
  estimated_entropy <- calculate_estimated_entropy(site_wise_info)
  
  if(give_p = TRUE){
    return( t.test(site_wise_info, mu = null_entropy, alternative = "less") )
  }  else{
  return(c(null_entropy, estimated_entropy))
  }
}
