################################################################################
################################################################################
#CODE USED TO COMPARE CPU TIME OF CLASSICAL AND PROJECTED PARAMETRIC BOOTSTRAP
#IN SECTION 7.2
################################################################################
################################################################################

rm(list=ls())

################################################################################
#CODE FOR TIMING THE CLASSICAL BOOTSTRAP
################################################################################
############################
##START TIMING HERE
ptm <- proc.time()
#############################
library("truncdist")
set.seed(12345)
B<-1000 #<_In the paper B=100,000
L<-0
U<-1
K<-100
const<-5
delta<-(U-L)/K
centers<-seq(L+delta/2,U-delta/2,by=delta) 
beta.true=0.5
st.dev=0.2
theta_true<-c(const,beta.true)
mx<-function(x)(U-L)*const*(dtrunc(x,a=L,b=U,spec="norm", mean=beta.true,sd=st.dev))
lambda<-function(x)dtrunc(x,a=L,b=U,spec="norm", mean=beta.true,sd=st.dev)
m<-function(theta)(U-L)*theta[1]*(dtrunc(centers,a=L,b=U,spec="norm", mean=theta[2],sd=st.dev))
mk_true<-m(theta_true)
rPOI<-function(mui){rpois(1, lambda=mui)}
rPOI<-Vectorize(rPOI,"mui")
nuk<-rPOI(mk_true)
g_parallel_Pearson<-function(observed,expected){(observed-expected)/expected}
int1<-function(theta1)integrate(function(x)dtrunc(x,a=L,b=U,spec="norm",mean=theta1,sd=st.dev)*(x-theta1),lower=L,upper=U)$value
mdot1_over_m<-function(theta1)((centers-theta1)-int1(theta1))/(st.dev^2)
est_eqMLE<-function(theta1,nu)sum(mdot1_over_m(theta1)*(nu-m(c(const,theta1))))
est_eqMLE<-Vectorize(est_eqMLE,"theta1")
beta_hat_obs<-uniroot(est_eqMLE, lower = 0.00001,upper=5,nu=nuk)$root
c_hat_obs<-mean(nuk)
theta_hat_obs<-c(c_hat_obs,beta_hat_obs)
mk_hat_obs<-m(theta_hat_obs)
KS_Linearized_hat<-c()
for(b in 1:B){
  nu0<-rPOI(mk_hat_obs)
  beta_hat<-uniroot(est_eqMLE, lower = 0.00001,upper=5,nu=nu0)$root
  c_hat<-mean(nu0)
  theta_hat<-c(c_hat,beta_hat)
  mk_hat<-m(theta_hat)
  vk_Linearized_hat<-cumsum(g_parallel_Pearson(nu0,mk_hat))/sqrt(K)
  KS_Linearized_hat[b]<-max(abs(vk_Linearized_hat))
  print(b)
}
############################
##END TIMING HERE
time1 <- proc.time()-ptm 
time1
#############################

tot1<-time1[1]+time1[2]
#> tot1
#user.self 
#18.38


rm(list=ls())

################################################################################
#CODE FOR TIMING THE PROJECTED BOOTSTRAP
################################################################################

############################
##START TIMING HERE
ptm <- proc.time()
#############################
library("truncdist")
set.seed(12345)
B<-1000 #<_In the paper B=100,000
L<-0
U<-1
K<-100
const<-5
delta<-(U-L)/K 
centers<-seq(L+delta/2,U-delta/2,by=delta)  
beta.true=0.5
st.dev=0.2
theta_true<-c(const,beta.true)
lambda<-function(x)dtrunc(x,a=L,b=U,spec="norm", mean=beta.true,sd=st.dev)
m<-function(theta)(U-L)*theta[1]*(dtrunc(centers,a=L,b=U,spec="norm", mean=theta[2],sd=st.dev))
mk_true<-m(theta_true)
rPOI<-function(mui){rpois(1, lambda=mui)}
rPOI<-Vectorize(rPOI,"mui")
nuk<-rPOI(mk_true)
g_parallel_Pearson<-function(observed,expected){(observed-expected)/expected}
int1<-function(theta1)integrate(function(x)dtrunc(x,a=L,b=U,spec="norm",mean=theta1,sd=st.dev)*(x-theta1),lower=L,upper=U)$value
mdot1_over_m<-function(theta1)((centers-theta1)-int1(theta1))/(st.dev^2)
est_eqMLE<-function(theta1,nu)sum(mdot1_over_m(theta1)*(nu-m(c(const,theta1))))
est_eqMLE<-Vectorize(est_eqMLE,"theta1")
beta_hat_obs<-uniroot(est_eqMLE, lower = 0.00001,upper=5,nu=nuk)$root
c_hat_obs<-mean(nuk)
theta_hat_obs<-c(c_hat_obs,beta_hat_obs)
mk_hat_obs<-m(theta_hat_obs)
mdot1_over_m_obs<-mdot1_over_m(beta_hat_obs)
coefficient1<-(cumsum(rep(1,length(centers))/c_hat_obs)/K)/mean(mk_hat_obs/(c_hat_obs^2))
coefficient2<-(cumsum(mdot1_over_m_obs)/K)/mean(mk_hat_obs*mdot1_over_m_obs^2)
KS_Linearized_tilde<-c()
for(b in 1:B){
  nu0<-rPOI(mk_hat_obs)
  vk_Linearized_tilde<-cumsum(g_parallel_Pearson(nu0,mk_hat_obs))/sqrt(K)-sqrt(K)*coefficient1*mean((nu0-mk_hat_obs)/c_hat_obs)-
    sqrt(K)*coefficient2*mean(mdot1_over_m_obs*(nu0-mk_hat_obs))
  KS_Linearized_tilde[b]<-max(abs(vk_Linearized_tilde))
  print(b)
}
############################
##END TIMING HERE
time2 <- proc.time()-ptm 
time2
#############################

tot2<-time2[1]+time2[2]
tot2
