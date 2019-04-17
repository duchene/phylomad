`tripleProb` <-
function (para)
{
   prob<-.C("tripleProb",as.double(para[1]),as.double(para[2]),as.double(para[3]),as.double(para[4]),p0=double(1),p1=double(1),p2=double(1),p3=double(1),p4=double(1),PACKAGE="phybase")

   p<-rep(-1,5)
   p[1] <- prob$p0
   p[2] <- prob$p1
   p[3] <- prob$p2
   p[4] <- prob$p3
   p[5] <- prob$p4
   return(p)
}