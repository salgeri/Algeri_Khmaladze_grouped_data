################################################################################
################################################################################
                      #Figure 3 (left panel) 
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

# g for Pearson linearized counterpart
g_parallel_Pearson<-function(observed,expected){(observed-expected)/expected}


#Estimating equations
int1<-function(theta1)integrate(function(x)dtrunc(x,a=L,b=U,spec="norm",mean=theta1,sd=st.dev)*(x-theta1),lower=L,upper=U)$value
mdot1_over_m<-function(theta1)((centers-theta1)-int1(theta1))/(st.dev^2)
est_eqMLE<-function(theta1,nu)sum(mdot1_over_m(theta1)*(nu-m(c(const,theta1))))
est_eqMLE<-Vectorize(est_eqMLE,"theta1")

#``Observed'' parameter values
beta_hat_obs<-uniroot(est_eqMLE, lower = 0.00001,upper=5,nu=nuk)$root
c_hat_obs<-mean(nuk)
theta_hat_obs<-c(c_hat_obs,beta_hat_obs)

#Quantities needed for evaluations of projected process
mk_hat_obs<-m(theta_hat_obs)
mdot1_over_m_obs<-mdot1_over_m(beta_hat_obs)
coefficient1<-(cumsum(rep(1,length(centers))/c_hat_obs)/K)/mean(mk_hat_obs/(c_hat_obs^2))
coefficient2<-(cumsum(mdot1_over_m_obs)/K)/mean(mk_hat_obs*mdot1_over_m_obs^2)

##########SIMULATION OF NULL DISTRIBUTION with CLASSICAL AND PROJECTED
##########PARAMETRIC BOOTSTRAP

#Number of bootstrap replicates
B<-1000 #<-In the paper B=100,000

KS_Linearized_hat<-KS_Linearized_bar<-c()
for(b in 1:B){#b=1

    nu0<-rPOI(mk_hat_obs)
  beta_hat<-uniroot(est_eqMLE, lower = 0.00001,upper=5,nu=nu0)$root
  c_hat<-mean(nu0)
  theta_hat<-c(c_hat,beta_hat)
  mk_hat<-m(theta_hat)
  
  #Estimated process
  vk_Linearized_hat<-cumsum(g_parallel_Pearson(nu0,mk_hat))/sqrt(K)
  KS_Linearized_hat[b]<-max(abs(vk_Linearized_hat))
  
  #Projected process
  vk_Linearized_bar<-cumsum(g_parallel_Pearson(nu0,mk_hat_obs))/sqrt(K)-sqrt(K)*coefficient1*mean((nu0-mk_hat_obs)/c_hat_obs)-
                                   sqrt(K)*coefficient2*mean(mdot1_over_m_obs*(nu0-mk_hat_obs))
  KS_Linearized_bar[b]<-max(abs(vk_Linearized_bar))
  print(b)
}

#Obtaining the null distributions
LL<-10000
cseq<-seq(0,10,length=LL)
p_hat<-c()
p_bar<-c()
for(i in 1:LL){
  p_hat[i]<-mean(KS_Linearized_hat<=cseq[i])
  p_bar[i]<-mean(KS_Linearized_bar<=cseq[i])
}

#####FIGURE 3 (left panel)
par(mfrow=c(1,1),mgp=c(2.2,0.7,0),mar=c(4,4.5,0.2,0.5),oma=c(0.1,0.1,0.1,0.1))
plot(cseq,p_hat,type="l",col="black",ylim=c(0,1),xlim=c(0,1.7),lty=1,cex.lab=1.5,cex.axis=1.5,
     main="", lwd=2,xlab=expression(s),ylab=expression(paste(P(KS<=s))))
lines(cseq,p_bar,type="l",lwd=2,lty=2,col="darkorange")
legend("bottomright",legend=c( expression(hat(KS)), expression(bar(KS))),
       lty=c( 1,2),#,3),
       lwd=c(rep(2,2)),#,4),
       col=c("black", "darkorange"),
       box.col =  rgb(0, 0, 0,0),cex=1.5)


