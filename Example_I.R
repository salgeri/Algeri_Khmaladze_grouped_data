################################################################################
################################################################################
##                                                                            ##
## Code needed to replicate the results obtained for Example I, i.e.,         ##
## the power simulation described in Sec 6.2 and left panel of Figure 2       ##
##                                                                            ##
################################################################################
################################################################################

#Packges used################################################
library("truncdist")
#############################################################

#Set the seed to ensure reproducibility
set.seed(123456)

#Region \X being considered
L<-0
U<-1

#Number of bins
K<-100

#constant to which T/K converges
const<-5
#Expected sample size
TT<-const*K

#We make the bins equally spaced over \X
delta<-(U-L)/K #Bin width
centers<-seq(L+delta/2,U-delta/2,by=delta)  #Centers of the bins

#######NULL MODEL

#true value of parameter beta
beta_true=1.5
#true value of parameter theta
theta_true<-c(const,beta_true)
#lambda density describing the spatial spread
lambda<-function(x)dtrunc(x,a=L,b=U,spec="exp", rate=beta_true)
#expected counts under the null as a function of theta
m<-function(theta)(U-L)*theta[1]*dtrunc(centers,a=L,b=U,spec="exp", rate=theta[2])
#expected counts under the null at the true value of theta
mk_true<-m(theta_true)

#Function generating observed Poisson counts given a vector of means
rPOI<-function(mui){rpois(1, lambda=mui)}
rPOI<-Vectorize(rPOI,"mui")

# g for Pearson chi-square
g_Pearson<-function(observed,expected){(observed-expected)^2/expected-1}

#######FUNCTIONAL DIRECTIONS DEFINING THE ALTERNATIVES

#######Defining h1(x)
rr<-integrate(function(x)log(x)*exp(-beta_true*x),lower=L,upper=U)$value/integrate(function(x)exp(-beta_true*x),lower=L,upper=U)$value
h0<-function(x)log(x)-rr
h02<-function(x)lambda(x)*(h0(x))^2
norm_h<-sqrt(integrate(h02,lower=L,upper=U)$value)
h1<-function(x)h0(x)/norm_h

#Constant ``a'' in the definition of h1 in the paper
1/norm_h
#Constant ``b'' in the definition of h1 in the paper
-rr/norm_h

#Alternative model
mk_H1<-mk_true*(1+h1(centers)/sqrt(K*const))

#######Defining hat(h1)(x)
#Score function for parametric family
phi_beta<-(1/beta_true-exp(-beta_true)/(1-exp(-beta_true)))
score<-function(x){phi_beta-x}
norm2_score<-integrate(function(x)(score(x))^2*lambda(x),lower=L,upper=U)$value
inner_score_int<-function(x)h1(x)*score(x)*lambda(x)
inner_score<-integrate(inner_score_int,lower=L,upper=U)$value
mean_h<-integrate(function(x)h1(x)*lambda(x),lower=L,upper=U)$value
h1_hat<-function(x)h1(x)-(inner_score/norm2_score)*score(x)-mean_h
shift<-integrate(h1_hat,lower=L,upper=U)$value/sqrt(const)
shift

#Estimating euqations
int1<-function(beta){integrate(function(x)x*exp(-beta*x),lower=L,upper=U)$value/integrate(function(x)exp(-beta*x),lower=L,upper=U)$value}
lambdadot1_over_lambda<-function(beta)(int1(beta)-centers)
est_eqMLE<-function(beta,chat,nu)sum(lambdadot1_over_lambda(beta)*(nu-m(c(chat,beta))))
est_eqMLE<-Vectorize(est_eqMLE,"beta")

#Number of bootstrap replicates
B<-1000 #<-In the paper B=100,000

#######POWER SIMULATION
Pearson_hat<-c()
Pearson1_hat<-c()

alpha=0.05 #Significance level

for(b in 1:B){
  
  #Simulation under $H_0$ to get the quantiles
  nu0<-rPOI(mk_true)
  c_hat<-sum(nu0)/K
  beta_hat<-uniroot(est_eqMLE, lower = 0.00001,upper=5,nu=nu0,chat=c_hat)$root
  theta_hat<-c(c_hat,beta_hat)
  mk_hat<-m(theta_hat)
  Pearson_hat[b]<-sum(g_Pearson(nu0,mk_hat))/sqrt(K)
  
  #Simulation under $H_1$ with $h_1$ to get the power
  nu1<-rPOI(mk_H1)
  c_hat1<-sum(nu1)/K
  beta_hat1<-uniroot(est_eqMLE, lower = 0.00001,upper=5,nu=nu1,chat=c_hat1)$root
  theta_hat1<-c(c_hat1,beta_hat1)
  mk_hat1<-m(theta_hat1)
  Pearson1_hat[b]<-sum(g_Pearson(nu1,mk_hat1))/sqrt(K)
  
  print(b)
}


sd(Pearson1_hat)

qL_Pearson_hat<-quantile(Pearson_hat,alpha/2)
qU_Pearson_hat<-quantile(Pearson_hat,1-alpha/2)

power_Pearson1_hat<-mean((Pearson1_hat<qL_Pearson_hat)|(Pearson1_hat>qU_Pearson_hat))

#SIMULATED POWER
power_Pearson1_hat

#FIGURE 2 LEFT PANEL
epsilon<-0.00001
par(mfrow=c(1,1),mgp=c(2.2,0.7,0),mar=c(4,3,0.2,0.5),oma=c(0.1,0.1,0.1,0.1))

m_tilde<-function(x)const*lambda(x)*(1+h1(x)/sqrt(TT))
m_tilde_hat<-function(x)const*lambda(x)*(1+h1_hat(x)/sqrt(TT))

plot(seq(L+epsilon,U-epsilon,by=0.01),const*lambda(seq(L+epsilon,U-epsilon,by=0.01)),type="l",col="black",
     lwd=3,ylim=c(2,10),
     cex.lab=1.5,cex.axis=1.5,
     xlab="x",ylab=" ")
lines(seq(L+epsilon,U-epsilon,by=0.01),m_tilde(seq(L+epsilon,U-epsilon,by=0.01)),col="dodgerblue1",lwd=3,lty=2)
lines(seq(L+epsilon,U-epsilon,by=0.01),m_tilde_hat(seq(L+epsilon,U-epsilon,by=0.01)),col="tomato3",lwd=3,lty=4)
legend("topright",legend=c(expression(m[theta](x)),expression(paste(m[theta](x),"[","1+",textstyle(frac(h[1](x),sqrt(T))),"]")), 
                           expression(paste(m[theta](x),"[1+",textstyle(frac(hat(h)[1](x),sqrt(T))),"]"))),
       lty=c( 1,2,4),
       lwd=c(rep(3,3)),#,4),
       col=c("black","dodgerblue1","tomato3"),
      box.col =  rgb(0, 0, 0,0),cex=1.5)
