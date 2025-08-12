################################################################################
################################################################################
#Derivation of the correction term 0.6/sqrt(K)

################################################################################
################################################################################

#Packges used################################################
library("truncdist")
#############################################################

#Set the seed to ensure reproducibility
set.seed(12345)

#Region \X being considered
L<-0
U<-1

#Maximum number of bins considered
K<-10000

#constant to which T/K converges
const<-5

#We make the bins equally spaced over \X
delta<-(U-L)/K #Bin width
centers<-seq(L+delta/2,U-delta/2,by=delta)  #Centers of the bins

#######NULL MODEL
mk_true<-rep(const,K) #We choose the men to be constant over all the bins

#Function generating observed Poisson counts given a vector of means
rPOI<-function(mui){rpois(1, lambda=mui)}
rPOI<-Vectorize(rPOI,"mui")

#Different number of bins considered
Kseq<-seq(10,K-10,by=10)

#Number of bootstrap replicates
B<-1000  #In the paper B=100,000 has been used 

#We now construct a process of partial sums which
#we know converges to a Brownian Bridge for large enough K
#(we simply chose the projected normalized linear statistics
#with mean equal to c over all the bins)
#and we proceed simulating the values of KS statistic 
#computed as the maximum of the absolute value of such a process
#for each value of K in Kseq

KS_ref<-matrix(rep(0,B*length(Kseq)),nrow=B,ncol=length(Kseq))
for(b in 1:B){#b=1
  nu0<-rPOI(mk_true)
  vk<-cumsum((nu0-mk_true)/sqrt(mk_true))/sqrt(K)-(cumsum(rep(1,K))/K)*sum((nu0-mk_true)/sqrt(mk_true))/sqrt(K)
  for(k in 1:length(Kseq)){
    vk_new<-vk[seq(1,K,length=Kseq[k])]
    KS_ref[b,k]<-max(abs(c(0,vk_new)))
  }
  print(b)
}

#Deriving the null distribution of the test statistic
#for each value of K considered
LL<-1000
cseq<-seq(0,3,length=LL)
p_ref<-matrix(rep(0,LL*length(Kseq)),nrow=LL,ncol=length(Kseq))
for(i in 1:LL){
  for(k in 1:length(Kseq)){
    p_ref[i,k]<-mean(KS_ref[,k]<=cseq[i])
  }}

#Kolmorogov distribution (theoretical limit)
my_pKS<-function(x){
  kvec<-1:100
  summands=exp(-((2*kvec-1)^2)*pi^2/(8*x^2))
  (sqrt(2*pi)/x)*sum(summands)
}
my_pKS<-Vectorize(my_pKS)

#Deriving the shift between the distribution
#of the KS statistic and the  Kolmogorov distribution
shifts<-seq(0.001,1,by=0.001)
diff_fun<-function(yy)apply(abs(p_ref-my_pKS(cseq+yy)),2,mean)
diff_fun<-Vectorize(diff_fun)
diff<-t(sapply(shifts,diff_fun))
shifts_values<-shifts[apply(diff,2,which.min)]

#Regressing the values of the shift over the values of K considered
betas<-lm(log(shifts_values)~log(Kseq))$coeff
exp(betas[1])
betas[2]

#> exp(betas[1])
#(Intercept) 
#0.56874325 
#> betas[2]
#log(Kseq) 
#-0.49719094 

#Graphical summary of the results
par(mfrow=c(1,1),mgp=c(2.2,0.7,0),mar=c(3.3,4.1,2,2),oma=c(0.1,0.1,0.1,0.1))
theoretical_shift<-function(K){0.6/sqrt(K)}
theoretical_shift<-Vectorize(theoretical_shift)
plot(Kseq,shifts_values, ylab="Shift",xlab="K",ylim=c(-0.001,0.2),xlim=c(10,10000),pch=16,lwd=2,
     cex.lab=1.5,cex.axis=1.5)
abline(h=0,col="grey67",lwd=2,lty=2)
legend("topright",legend=expression(0.6/sqrt(K)),
       lty=c(1),#,3),
       lwd=c(2),#,4),
       col=c("tomato3"),box.col =  rgb(0, 0, 0,0),cex=1.5)
lines(seq(9.2,10000,by=0.01),theoretical_shift(seq(9.2,10000,by=0.01)),col="tomato3",lwd=2,lty=1)

