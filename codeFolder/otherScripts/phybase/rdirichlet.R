'rdirichlet'<-
function(n,a) 
## pick n random deviates from the Dirichlet function with shape 
## parameters a 
{ 
    l<-length(a); 
    x<-matrix(rgamma(l*n,a),ncol=l,byrow=TRUE); 
    sm<-x%*%rep(1,l); 
    x/as.vector(sm); 
} 
