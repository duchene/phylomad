`simnucleotide` <-
function(nucleotide, brlens, qmatrix)
{
	tranprob<-expm(Matrix(qmatrix*brlens))
	new<-tranprob[nucleotide,]
	newnucleotide<-apply(new,1,function(x) sample(1:4,1,TRUE,x))
	newnucleotide
}

