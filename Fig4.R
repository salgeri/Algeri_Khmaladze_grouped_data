################################################################################
################################################################################
#Code needed to replicate Figure 4
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

#Number of bins
K<-50  #<------In the paper K=50,100,1000 were used

#constant to which T/K converges
const<-5

#We make the bins equally spaced over \X
delta<-(U-L)/K #Bin width
centers<-seq(L+delta/2,U-delta/2,by=delta)  #Centers of the bins

#NORMAL MODEL
#true value of the parameters
beta.true=0.5
st.dev=0.2
theta_true<-c(const,beta.true)
#lambda density describing the spatial spread
lambda<-function(x)dtrunc(x,a=L,b=U,spec="norm", mean=beta.true,sd=st.dev)
#expected counts under the null as a function of theta
m<-function(theta)(U-L)*theta[1]*(dtrunc(centers,a=L,b=U,spec="norm", mean=theta[2],sd=st.dev))
#expected counts under the null at the true value of theta
mk_true<-m(theta_true)

#EXPONENTIAL MODEL
#true value unknown parameters
beta.exp.true=0.5
theta.exp_true<-c(const,beta.exp.true)
#lambda density describing the spatial spread
lambda.exp<-function(x)dtrunc(x,a=L,b=U,spec="exp", rate=beta.exp.true)
#expected counts under the null as a function of theta
m.exp<-function(theta)(U-L)*theta[1]*(dtrunc(centers,a=L,b=U,spec="exp", rate=theta[2]))
#expected counts under the null at the true value of theta
mk.exp_true<-m.exp(theta.exp_true)

#Function generating observed Poisson counts given a vector of means
rPOI<-function(mui){rpois(1, lambda=mui)}
rPOI<-Vectorize(rPOI,"mui")

#``Observed'' counts 
nuk<-rPOI(mk_true)
nuk.exp<-rPOI(mk.exp_true)

#Estimating equations Normal model and estimates of the parameters
int1<-function(theta1)integrate(function(x)dtrunc(x,a=L,b=U,spec="norm",mean=theta1,sd=st.dev)*(x-theta1),lower=L,upper=U)$value
mdot1_over_m<-function(theta1)((centers-theta1)-int1(theta1))/(st.dev^2)
est_eqMLE<-function(theta1,nu)sum(mdot1_over_m(theta1)*(nu-m(c(const,theta1))))
est_eqMLE<-Vectorize(est_eqMLE,"theta1")
beta_hat_obs<-uniroot(est_eqMLE, lower = 0.00001,upper=5,nu=nuk)$root
c_hat_obs<-mean(nuk)
theta_hat_obs<-c(c_hat_obs,beta_hat_obs)
mk_hat_obs<-m(theta_hat_obs)
mdot1_over_m_obs<-mdot1_over_m(beta_hat_obs)

#Estimating equations Exponential model and estimates of the parameters
int2<-function(theta1)integrate(function(x)dtrunc(x,a=L,b=U,spec="exp",rate=theta1)*x,lower=L,upper=U)$value
mdot1_over_m.exp<-function(theta1)(int2(theta1)-centers)
est_eqMLE.exp<-function(theta1,nu)sum(mdot1_over_m.exp(theta1)*(nu-m.exp(c(const,theta1))))
est_eqMLE.exp<-Vectorize(est_eqMLE.exp,"theta1")
beta.exp_hat_obs<-uniroot(est_eqMLE.exp, lower = 0.00001,upper=10,nu=nuk.exp)$root
theta.exp_hat_obs<-c(c_hat_obs,beta.exp_hat_obs)
mk.exp_hat_obs<-m.exp(theta.exp_hat_obs)
mdot1_over_m.exp_obs<-mdot1_over_m.exp(beta.exp_hat_obs)

#Quantities needed for constructing the transformed process
#based on the U_p operator for both models
#To speed up the simualtion, the transformation 
#has been derived explicitly and all the inner products involved are 
#computed below

#Indicator functions defining the sets B1 and B2
I_B1<-as.numeric(centers<=0.5)
I_B2<-as.numeric(centers>0.5)

#Reference functions r_j, j=1,2, for both models
r1<-function(nu)I_B1*(nu-mk_hat_obs)/sqrt(mean(I_B1)*mk_hat_obs)
r2<-function(nu)I_B2*(nu-mk_hat_obs)/sqrt(mean(I_B2)*mk_hat_obs)
r1.exp<-function(nu)I_B1*(nu-mk.exp_hat_obs)/sqrt(mean(I_B1)*mk.exp_hat_obs)
r2.exp<-function(nu)I_B2*(nu-mk.exp_hat_obs)/sqrt(mean(I_B2)*mk.exp_hat_obs)

#Normalized Score function s_j, j=1,2, for both models
s1<-function(nu)(nu-mk_hat_obs)/sqrt(mean(mk_hat_obs))
s2<-function(nu)mdot1_over_m_obs*(nu-mk_hat_obs)/sqrt(mean(mk_hat_obs*mdot1_over_m_obs^2))
s1.exp<-function(nu)(nu-mk.exp_hat_obs)/sqrt(mean(mk.exp_hat_obs))
s2.exp<-function(nu)mdot1_over_m.exp_obs*(nu-mk.exp_hat_obs)/sqrt(mean(mk.exp_hat_obs*mdot1_over_m.exp_obs^2))

#summands of normalized linear statistic for both models
lg_parallel<-function(nu){(nu-mk_hat_obs)/sqrt(mk_hat_obs)}
lg_parallel.exp<-function(nu){(nu-mk.exp_hat_obs)/sqrt(mk.exp_hat_obs)}

#Inner products <\ell_t,r_j>, j=1,2
inner_lgA_r1<-inner_lgA_r1.exp<-(cumsum(I_B1)/K)/sqrt(mean(I_B1))
inner_lgA_r2<-inner_lgA_r2.exp<-(cumsum(I_B2)/K)/sqrt(mean(I_B2))

#Inner products <\ell_t,s_j>, j=1,2, for both models
inner_lgA_s1<-(cumsum(sqrt(mk_hat_obs))/K)/sqrt(mean(mk_hat_obs))
inner_lgA_s2<-(cumsum(mdot1_over_m_obs*sqrt(mk_hat_obs))/K)/sqrt(mean(mk_hat_obs*mdot1_over_m_obs^2))
inner_lgA_s1.exp<-(cumsum(sqrt(mk.exp_hat_obs))/K)/sqrt(mean(mk.exp_hat_obs))
inner_lgA_s2.exp<-(cumsum(mdot1_over_m.exp_obs*sqrt(mk.exp_hat_obs))/K)/sqrt(mean(mk.exp_hat_obs*mdot1_over_m.exp_obs^2))

#Inner products <r_j,s_1>, j=1,2, for both models (needed to construct tilde(r)_j functions)
inner_r1s1<-mean(I_B1*sqrt(mk_hat_obs))/sqrt(mean(I_B1)*mean(mk_hat_obs))
inner_r2s1<-mean(I_B2*sqrt(mk_hat_obs))/sqrt(mean(I_B2)*mean(mk_hat_obs))
inner_r1s1.exp<-mean(I_B1*sqrt(mk.exp_hat_obs))/sqrt(mean(I_B1)*mean(mk.exp_hat_obs))
inner_r2s1.exp<-mean(I_B2*sqrt(mk.exp_hat_obs))/sqrt(mean(I_B2)*mean(mk.exp_hat_obs))

#Inner products <r_j,s_2>, j=1,2, for both models (needed to construct tilde(r)_j functions)
inner_r1s2<-mean(I_B1*sqrt(mk_hat_obs)*mdot1_over_m_obs)/sqrt(mean(I_B1)*mean(mk_hat_obs*mdot1_over_m_obs^2))
inner_r2s2<-mean(I_B2*sqrt(mk_hat_obs)*mdot1_over_m_obs)/sqrt(mean(I_B2)*mean(mk_hat_obs*mdot1_over_m_obs^2))
inner_r1s2.exp<-mean(I_B1*sqrt(mk.exp_hat_obs)*mdot1_over_m.exp_obs)/sqrt(mean(I_B1)*mean(mk.exp_hat_obs*mdot1_over_m.exp_obs^2))
inner_r2s2.exp<-mean(I_B2*sqrt(mk.exp_hat_obs)*mdot1_over_m.exp_obs)/sqrt(mean(I_B2)*mean(mk.exp_hat_obs*mdot1_over_m.exp_obs^2))

#Constructing tilde(r)_j functions for both models
r2_tilde<-function(nu)r2(nu)-inner_r2s1*(s1(nu)-r1(nu))/(1-inner_r1s1)
r2_tilde.exp<-function(nu)r2.exp(nu)-inner_r2s1.exp*(s1.exp(nu)-r1.exp(nu))/(1-inner_r1s1.exp)

#Inner products <\ell_t,tilde(r)_2> for both models
inner_lgA_r2_tilde<-inner_lgA_r2-(inner_r2s1)*(inner_lgA_s1-inner_lgA_r1)/(1-inner_r1s1)
inner_lgA_r2_tilde.exp<-inner_lgA_r2.exp-(inner_r2s1.exp)*(inner_lgA_s1.exp-inner_lgA_r1.exp)/(1-inner_r1s1.exp)

#Inner products <s_2,tilde(r)_2> for both models
inner_s2_r2_tilde<-inner_r2s2+inner_r2s1*inner_r1s2/(1-inner_r1s1)
inner_s2_r2_tilde.exp<-inner_r2s2.exp+inner_r2s1.exp*inner_r1s2.exp/(1-inner_r1s1.exp)

#Terms involved in the U_p transformation and based on the inner products above
W1<-(inner_lgA_s1-inner_lgA_r1)/(1-inner_r1s1)
W2<-(inner_lgA_s2-inner_lgA_r2_tilde)/(1-inner_s2_r2_tilde)
Z12<-(inner_r2s1-inner_r1s2)/(1-inner_s2_r2_tilde)
W1.exp<-(inner_lgA_s1.exp-inner_lgA_r1.exp)/(1-inner_r1s1.exp)
W2.exp<-(inner_lgA_s2.exp-inner_lgA_r2_tilde.exp)/(1-inner_s2_r2_tilde.exp)
Z12.exp<-(inner_r2s1.exp-inner_r1s2.exp)/(1-inner_s2_r2_tilde.exp)


#Number of bootstrap replicates
B<-1000  #In the paper B=100,000 has been used 

#Simulation of the distribution of the bar(KS) and the KS* statistics
#under the null for the Exponential and normal model
mat_U<-mat_U.exp<-matrix(rep(0,B*K),nrow=B,ncol=K)
KS<-KS.exp<-KS_U<-KS_U.exp<-c()
for(b in 1:B){
  
  ##############NORMAL MODEL
  nu0<-rPOI(mk_hat_obs)
  
  #Transformed process and KS*
  xiA<-cumsum(lg_parallel(nu0))/sqrt(K)-
    inner_lgA_r1*sqrt(K)*mean(s1(nu0))-
    inner_lgA_r2*sqrt(K)*mean(s2(nu0))
  vk_U<-xiA-W1*sqrt(K)*mean((s1(nu0)-r1(nu0)))-(W2-W1*Z12)*sqrt(K)*mean((s2(nu0)-r2_tilde(nu0)))
  mat_U[b,]<-vk_U
  KS_U[b]<-max(abs(vk_U))
  
  #Projected process and bar(KS)
  vk<-cumsum(lg_parallel(nu0))/sqrt(K)-
    inner_lgA_s1*sqrt(K)*mean(s1(nu0))-
    inner_lgA_s2*sqrt(K)*mean(s2(nu0))
  KS[b]<-max(abs(vk))
  
  
  ##############EPONENTIAL MODEL
  nu0.exp<-rPOI(mk.exp_hat_obs)
  
  #Transformed process and KS*
  xiA.exp<-cumsum(lg_parallel.exp(nu0.exp))/sqrt(K)-
    inner_lgA_r1.exp*sqrt(K)*mean(s1.exp(nu0.exp))-
    inner_lgA_r2.exp*sqrt(K)*mean(s2.exp(nu0.exp))
  vk_U.exp<-xiA.exp-W1.exp*sqrt(K)*mean((s1.exp(nu0.exp)-r1.exp(nu0.exp)))-
    (W2.exp-W1.exp*Z12.exp)*sqrt(K)*mean((s2.exp(nu0.exp)-r2_tilde.exp(nu0.exp)))
  mat_U.exp[b,]<-vk_U.exp
  KS_U.exp[b]<-max(abs(vk_U.exp))
  
  #Projected process and bar(KS)
  vk.exp<-cumsum(lg_parallel.exp(nu0.exp))/sqrt(K)-
    inner_lgA_s1.exp*sqrt(K)*mean(s1.exp(nu0.exp))-
    inner_lgA_s2.exp*sqrt(K)*mean(s2.exp(nu0.exp))
  KS.exp[b]<-max(abs(vk.exp))
  
  print(c(b,KS_U[b],KS_U.exp[b]))
}


#Deriving the null distribution of the bar(KS) and KS* statistics
LL<-1000
cseq<-seq(0,3,length=LL)
pp<-pp.exp<-c()
p_U<-p_U.exp<-c()
for(i in 1:LL){
  pp[i]<-mean(KS<=cseq[i])
  pp.exp[i]<-mean(KS.exp<=cseq[i])
  p_U[i]<-mean(KS_U+0.6/sqrt(K)<=cseq[i])
  p_U.exp[i]<-mean(KS_U.exp+0.6/sqrt(K)<=cseq[i])
  }

#Kolmorogov distribution
my_pKS<-function(x){
  kvec<-1:100
  summands=exp(-((2*kvec-1)^2)*pi^2/(8*x^2))
  (sqrt(2*pi)/x)*sum(summands)
}
my_pKS<-Vectorize(my_pKS)


########FIGURE 4 top left panel
par(mfrow=c(1,1),mgp=c(2.2,0.7,0),mar=c(3.3,4.1,2,0.5),oma=c(0.1,0.1,0.1,0.1))
plot(cseq,(my_pKS(cseq*sqrt(2)))^2,type="l",col="black",ylim=c(0,1),xlim=c(0.25,1.3),lty=1,cex.lab=1.5,cex.axis=1.5,cex.main=1.5,
     main="K=50", lwd=2,xlab=expression(s),ylab=expression(paste(P(bar(KS)<=s))))
lines(cseq,pp.exp,type="l",lwd=2,lty=4,col="dodgerblue")
lines(cseq,pp,type="l",lwd=2,lty=2,col="darkorange")

legend("bottomright",legend=c("Limit", "Exponential", "Normal"),
       lty=c(1,2,4),#,3),
       lwd=c(rep(2,3)),#,4),
       col=c("black", "darkorange","dodgerblue"),box.col =  rgb(0, 0, 0,0),cex=1.5)



########FIGURE 4 bottom left panel
par(mfrow=c(1,1),mgp=c(2.2,0.7,0),mar=c(3.3,4.1,2,0.5),oma=c(0.1,0.1,0.1,0.1))
plot(cseq,(my_pKS(cseq*sqrt(2)))^2,type="l",col="black",ylim=c(0,1),xlim=c(0.25,1.3),lty=1,cex.lab=1.5,cex.axis=1.5,cex.main=1.5,
     main=paste("K=",K), lwd=2,xlab=expression(s),ylab=expression(paste(P(KS^"*"+0.6/sqrt(K)<=s))))
lines(cseq,p_U.exp,type="l",lwd=2,lty=4,col="dodgerblue")
lines(cseq,p_U,type="l",lwd=2,lty=2,col="darkorange")

legend("bottomright",legend=c("Limit", "Exponential", "Normal"),
       lty=c(1,2,4),#,3),
       lwd=c(rep(2,3)),#,4),
       col=c("black", "darkorange","dodgerblue"),box.col =  rgb(0, 0, 0,0),cex=1.5)

