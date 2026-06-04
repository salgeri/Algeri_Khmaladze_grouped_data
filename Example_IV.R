################################################################################
################################################################################
        #Power simulation for Example IV  
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
K<-100

#constant to which T/K converges
const<-5

#We make the bins equally spaced over \X
delta<-(U-L)/K #Bin width
centers<-seq(L+delta/2,U-delta/2,by=delta)  #Centers of the bins


#######NULL MODEL

#true value parameters
beta_true=0.5
st.dev=0.2
theta_true<-c(const,beta_true)
#lambda density describing the spatial spread
lambda<-function(x)dtrunc(x,a=L,b=U,spec="norm", mean=beta_true,sd=st.dev)
#null_model as a function of x
mx<-function(x)const*(U-L)*lambda(x)
#expected counts under the null as a function of theta
m<-function(theta)(U-L)*theta[1]*(dtrunc(centers,a=L,b=U,spec="norm", mean=theta[2],sd=st.dev))
#expected counts under the null at the true value of theta
mk_true<-m(theta_true)

#Function generating observed Poisson counts given a vector of means
rPOI<-function(mui){rpois(1, lambda=mui)}
rPOI<-Vectorize(rPOI,"mui")
nuk<-rPOI(mk_true)

#Functions computing the Poisson cdf and pmf at a value r given a vector of means
#(to be used for the spectrum)
r=1  
pPOI<-function(mui){ppois(r, lambda=mui)}
pPOI<-Vectorize(pPOI,"mui")

dPOI<-function(mui){dpois(r-1, lambda=mui)}
dPOI<-Vectorize(dPOI,"mui")

# g for Pearson chi-square and the corresponding weighted linear statistic
g_Pearson<-function(observed,expected){(observed-expected)^2/expected-1}
g_parallel_Pearson<-function(observed,expected){(observed-expected)/expected}

g_spectrum<-function(observed,expected){ifelse(observed<=r,1,0)-pPOI(expected)}
g_parallel_spectrum<-function(observed,expected){dPOI(expected)*(observed-expected)}

#######FUNCTIONAL DIRECTIONS DEFINING THE ALTERNATIVES

#######Defining h4(x)
int2<-integrate(function(x)lambda(x)*(x-beta_true)^2,lower=L,upper=U)$value
h4<-function(x)((x-beta_true)^2-int2)/(2*st.dev^4)
norm_h4<-sqrt(integrate(function(x)lambda(x)*h4(x)^2,lower=L,upper=U)$value)
h4_norm<-function(x)h4(x)*(1/norm_h4)

#Constant ``a'' in the definition of h4 in the paper
(1/norm_h4)
#Constant ``b'' in the definition of h4 in the paper
int2

#alternative model as a function of x
mx_H1<-function(x)mx(x)*(1+h4_norm(x)/sqrt(const*K))
#expected counts under alternative on each bin
mk_H1<-mk_true*(1+h4_norm(centers)/sqrt(K*const))

#Comparing null and alternative
par(mfrow=c(1,1),mgp=c(2.2,0.7,0),mar=c(4,3,0.2,0.5),oma=c(0.1,0.1,0.1,0.1))

plot(centers,mk_true,type="l",col="black",
     lwd=2,ylim=c(0,10.5),
     cex.lab=1.5,cex.axis=1.5,
     xlab="x",ylab=" ")
lines(centers,mk_H1,col="tomato3",lwd=2,lty=2)
legend("topleft",legend=c(expression(m[theta](x)), expression(tilde(m)[theta](x))),
       lty=c( 1,2),
       lwd=c(rep(2,2)),#,4),
       col=c("black", "tomato3"),
       box.col =  rgb(0, 0, 0,0),cex=1.5)

#Estimating equations
int1<-function(theta1)integrate(function(x)dtrunc(x,a=L,b=U,spec="norm",mean=theta1,sd=st.dev)*(x-theta1),lower=L,upper=U)$value
int1(beta_true)
mdot1_over_m<-function(theta1)((centers-theta1)-int1(theta1))/(st.dev^2)
est_eqMLE<-function(theta1,nu,chat)sum(mdot1_over_m(theta1)*(nu-m(c(chat,theta1))))
est_eqMLE<-Vectorize(est_eqMLE,"theta1")

#Number of bootstrap replicates
B<-1000 #<-In the paper B=100,000
#######POWER SIMULATION
Pearson<-Linearized<-KS_Pearson<-KS_Linearized<-c()
Pearson1<-Linearized1<-KS_Pearson1<-KS_Linearized1<-c()
Pearson_hat<-Linearized_hat<-KS_Pearson_hat<-KS_Linearized_hat<-c()
Pearson1_hat<-Linearized1_hat<-KS_Pearson1_hat<-KS_Linearized1_hat<-c()

spectrum<-spectrum_Linearized<-KS_spectrum<-KS_spectrum_Linearized<-c()
spectrum1<-spectrum_Linearized1<-KS_spectrum1<-KS_spectrum_Linearized1<-c()
spectrum_hat<-spectrum_Linearized_hat<-KS_spectrum_hat<-KS_spectrum_Linearized_hat<-c()
spectrum1_hat<-spectrum_Linearized1_hat<-KS_spectrum1_hat<-KS_spectrum_Linearized1_hat<-c()
alpha=0.05
for(b in 1:B){#b=1
  #Simulation under $H_0$ to get the quantiles
  nu0<-rPOI(mk_true)
  Pearson[b]<-sum(g_Pearson(nu0,mk_true))/sqrt(K)
  Linearized[b]<-sum(g_parallel_Pearson(nu0,mk_true))/sqrt(K)
  vk_Pearson<-cumsum(g_Pearson(nu0,mk_true))/sqrt(K)
  vk_Linearized<-cumsum(g_parallel_Pearson(nu0,mk_true))/sqrt(K)
  KS_Pearson[b]<-max(abs(vk_Pearson))
  KS_Linearized[b]<-max(abs(vk_Linearized))
  spectrum[b]<-sum(g_spectrum(nu0,mk_true))/sqrt(K)
  spectrum_Linearized[b]<-sum(g_parallel_spectrum(nu0,mk_true))/sqrt(K)
  vk_spectrum<-cumsum(g_spectrum(nu0,mk_true))/sqrt(K)
  vk_spectrum_Linearized<-cumsum(g_parallel_spectrum(nu0,mk_true))/sqrt(K)
  KS_spectrum[b]<-max(abs(vk_spectrum))
  KS_spectrum_Linearized[b]<-max(abs(vk_spectrum_Linearized))
  
  #Estimating parameters
  c_hat<-sum(nu0)/K
  beta_hat<-uniroot(est_eqMLE, lower = 0.00001,upper=5,nu=nu0,chat=c_hat)$root
  theta_hat<-c(c_hat,beta_hat)
  mk_hat<-m(theta_hat)
  Pearson_hat[b]<-sum(g_Pearson(nu0,mk_hat))/sqrt(K)
  Linearized_hat[b]<-sum(g_parallel_Pearson(nu0,mk_hat))/sqrt(K)
  vk_Pearson_hat<-cumsum(g_Pearson(nu0,mk_hat))/sqrt(K)
  vk_Linearized_hat<-cumsum(g_parallel_Pearson(nu0,mk_hat))/sqrt(K)
  KS_Pearson_hat[b]<-max(abs(vk_Pearson_hat))
  KS_Linearized_hat[b]<-max(abs(vk_Linearized_hat))
  spectrum_hat[b]<-sum(g_spectrum(nu0,mk_hat))/sqrt(K)
  spectrum_Linearized_hat[b]<-sum(g_parallel_spectrum(nu0,mk_hat))/sqrt(K)
  vk_spectrum_hat<-cumsum(g_spectrum(nu0,mk_hat))/sqrt(K)
  vk_spectrum_Linearized_hat<-cumsum(g_parallel_spectrum(nu0,mk_hat))/sqrt(K)
  KS_spectrum_hat[b]<-max(abs(vk_spectrum_hat))
  KS_spectrum_Linearized_hat[b]<-max(abs(vk_spectrum_Linearized_hat))
  
  
  #Simulation under $H_1$ with $h_4$ to get the power
  nu1<-rPOI(mk_H1)
  Pearson1[b]<-sum(g_Pearson(nu1,mk_true))/sqrt(K)
  Linearized1[b]<-sum(g_parallel_Pearson(nu1,mk_true))/sqrt(K)
  vk_Pearson1<-cumsum(g_Pearson(nu1,mk_true))/sqrt(K)
  vk_Linearized1<-cumsum(g_parallel_Pearson(nu1,mk_true))/sqrt(K)
  KS_Pearson1[b]<-max(abs(vk_Pearson1))
  KS_Linearized1[b]<-max(abs(vk_Linearized1))
  spectrum1[b]<-sum(g_spectrum(nu1,mk_true))/sqrt(K)
  spectrum_Linearized1[b]<-sum(g_parallel_spectrum(nu1,mk_true))/sqrt(K)
  vk_spectrum1<-cumsum(g_spectrum(nu1,mk_true))/sqrt(K)
  vk_spectrum_Linearized1<-cumsum(g_parallel_spectrum(nu1,mk_true))/sqrt(K)
  KS_spectrum1[b]<-max(abs(vk_spectrum1))
  KS_spectrum_Linearized1[b]<-max(abs(vk_spectrum_Linearized1))
  
  #Estimating parameters
  c_hat1<-sum(nu1)/K
  beta_hat1<-uniroot(est_eqMLE, lower = 0.00001,upper=5,nu=nu1,chat=c_hat1)$root
  theta_hat1<-c(c_hat1,beta_hat1)
  mk_hat1<-m(theta_hat1)
  Pearson1_hat[b]<-sum(g_Pearson(nu1,mk_hat1))/sqrt(K)
  Linearized1_hat[b]<-sum(g_parallel_Pearson(nu1,mk_hat1))/sqrt(K)
  vk_Pearson1_hat<-cumsum(g_Pearson(nu1,mk_hat1))/sqrt(K)
  vk_Linearized1_hat<-cumsum(g_parallel_Pearson(nu1,mk_hat1))/sqrt(K)
  KS_Pearson1_hat[b]<-max(abs(vk_Pearson1_hat))
  KS_Linearized1_hat[b]<-max(abs(vk_Linearized1_hat))
  spectrum1_hat[b]<-sum(g_spectrum(nu1,mk_hat1))/sqrt(K)
  spectrum_Linearized1_hat[b]<-sum(g_parallel_spectrum(nu1,mk_hat1))/sqrt(K)
  vk_spectrum1_hat<-cumsum(g_spectrum(nu1,mk_hat1))/sqrt(K)
  vk_spectrum_Linearized1_hat<-cumsum(g_parallel_spectrum(nu1,mk_hat1))/sqrt(K)
  KS_spectrum1_hat[b]<-max(abs(vk_spectrum1_hat))
  KS_spectrum_Linearized1_hat[b]<-max(abs(vk_spectrum_Linearized1_hat))
  print(b)
}
#alpha<-0.1
qL_Pearson<-quantile(Pearson,alpha/2)
qU_Pearson<-quantile(Pearson,1-alpha/2)
qL_Linearized<-quantile(Linearized,alpha/2)
qU_Linearized<-quantile(Linearized,1-alpha/2)
qL_spectrum<-quantile(spectrum,alpha/2)
qU_spectrum<-quantile(spectrum,1-alpha/2)
qL_spectrum_Linearized<-quantile(spectrum_Linearized,alpha/2)
qU_spectrum_Linearized<-quantile(spectrum_Linearized,1-alpha/2)
q_KS_Pearson<-quantile(KS_Pearson,1-alpha)
q_KS_Linearized<-quantile(KS_Linearized,1-alpha)
q_KS_spectrum<-quantile(KS_spectrum,1-alpha)
q_KS_spectrum_Linearized<-quantile(KS_spectrum_Linearized,1-alpha)


qL_Pearson_hat<-quantile(Pearson_hat,alpha/2)
qU_Pearson_hat<-quantile(Pearson_hat,1-alpha/2)
qL_Linearized_hat<-quantile(Linearized_hat,alpha/2)
qU_Linearized_hat<-quantile(Linearized_hat,1-alpha/2)
qL_spectrum_hat<-quantile(spectrum_hat,alpha/2)
qU_spectrum_hat<-quantile(spectrum_hat,1-alpha/2)
qL_spectrum_Linearized_hat<-quantile(spectrum_Linearized_hat,alpha/2)
qU_spectrum_Linearized_hat<-quantile(spectrum_Linearized_hat,1-alpha/2)
q_KS_Pearson_hat<-quantile(KS_Pearson_hat,1-alpha)
q_KS_Linearized_hat<-quantile(KS_Linearized_hat,1-alpha)
q_KS_spectrum_hat<-quantile(KS_spectrum_hat,1-alpha)
q_KS_spectrum_Linearized_hat<-quantile(KS_spectrum_Linearized_hat,1-alpha)

power_Pearson1<-mean((Pearson1<qL_Pearson)|(Pearson1>qU_Pearson))
power_Linearized1<-mean((Linearized1<qL_Linearized)|(Linearized1>qU_Linearized))
power_spectrum1<-mean((spectrum1<qL_spectrum)|(spectrum1>qU_spectrum))
power_spectrum_Linearized1<-mean((spectrum_Linearized1<qL_spectrum_Linearized)|(spectrum_Linearized1>qU_spectrum_Linearized))
power_KS_Pearson1<-mean(KS_Pearson1>q_KS_Pearson)
power_KS_Linearized1<-mean(KS_Linearized1>q_KS_Linearized)
power_KS_spectrum1<-mean(KS_spectrum1>q_KS_spectrum)
power_KS_spectrum_Linearized1<-mean(KS_spectrum_Linearized1>q_KS_spectrum_Linearized)

power_Pearson1_hat<-mean((Pearson1_hat<qL_Pearson_hat)|(Pearson1_hat>qU_Pearson_hat))
power_Linearized1_hat<-mean((Linearized1_hat<qL_Linearized_hat)|(Linearized1_hat>qU_Linearized_hat))
power_spectrum1_hat<-mean((spectrum1_hat<qL_spectrum_hat)|(spectrum1_hat>qU_spectrum_hat))
power_spectrum_Linearized1_hat<-mean((spectrum_Linearized1_hat<qL_spectrum_Linearized_hat)|(spectrum_Linearized1_hat>qU_spectrum_Linearized_hat))
power_KS_Pearson1_hat<-mean(KS_Pearson1_hat>q_KS_Pearson_hat)
power_KS_Linearized1_hat<-mean(KS_Linearized1_hat>q_KS_Linearized_hat)
power_KS_spectrum1_hat<-mean(KS_spectrum1_hat>q_KS_spectrum_hat)
power_KS_spectrum_Linearized1_hat<-mean(KS_spectrum_Linearized1_hat>q_KS_spectrum_Linearized_hat)

#Power when testing simple hypotheses
power_Pearson1
power_Linearized1
power_spectrum1
power_spectrum_Linearized1
power_KS_Pearson1
power_KS_Linearized1
power_KS_spectrum1
power_KS_spectrum_Linearized1

#Power when testing parametric hypotheses
power_Pearson1_hat
power_Linearized1_hat
power_spectrum1_hat
power_spectrum_Linearized1_hat
power_KS_Pearson1_hat
power_KS_Linearized1_hat
power_KS_spectrum1_hat
power_KS_spectrum_Linearized1_hat


#When B=100,000
#> #Power when testing simple hypotheses
#  > power_Pearson1
#[1] 0.07161
#> power_Linearized1
#[1] 0.1053
#> power_spectrum1
#[1] 0.07796
#> power_spectrum_Linearized1
#[1] 0.12528
#> power_KS_Pearson1
#[1] 0.09056
#> power_KS_Linearized1
#[1] 0.11717
#> power_KS_spectrum1
#[1] 0.09578
#> power_KS_spectrum_Linearized1
#[1] 0.15216
#> 
#Power when testing parametric hypotheses
# power_Pearson1_hat
#[1] 0.07382
# power_Linearized1_hat
#[1] 0.1409
#> power_spectrum1_hat
#[1] 0.10162
#> power_spectrum_Linearized1_hat
#[1] 0.13492
#> power_KS_Pearson1_hat
#[1] 0.08283
#> power_KS_Linearized1_hat
#[1] 0.16466
#> power_KS_spectrum1_hat
#[1] 0.10096
#> power_KS_spectrum_Linearized1_hat
#[1] 0.15443
