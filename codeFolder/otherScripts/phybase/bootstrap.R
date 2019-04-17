`bootstrap` <-
function(sequence)
{ 
    nchar<-dim(sequence)[2]
    a<-sample(1:nchar,nchar,replace=TRUE)
    boot<-sequence[,a]
    boot
}