################################################################################
################################################################################
##                                                                            ##
## Code needed to reproduce the analysis of the X-ray source-free spectrum    ##
## from the Chandra X-ray observatory                                         ##
##                                                                            ##
################################################################################
################################################################################

#######################Installing and loading the required packages#############
#install.packages("rootSolve")
require(rootSolve)
################################################################################

#Setting the side for the numerical simulations
set.seed(1234)

#Loading data 


bkg<-read.table("C:\\Codes_Algeri_Khm_JRSSB2025\\Chandra_data.txt")

#Lower and upper bound of the interval defining the spectrum considered
L=14.6;U=17.4

#Values x_k corresponding to the midpoins of the bins
x<-bkg$bin_midpoints

#Bin frequencies
nu_obs<-bkg$counts

#Number of bins
K<-length(x)

#Fitted constant mean function
c_hat_obs<-sum(nu_obs)/K

#Mean function consisting of a linear model
m<-function(theta){theta[1]+theta[2]*x}

#Mean function consisting of a piecewise linear modelwith breakpoint at 15.6
break_point<-15.6
indicator<-ifelse(x<=break_point,1,0)
m2<-function(theta){theta[1]+theta[2]*(x-break_point)*indicator}

#Defining the estimating equations for theta in both models
est_eqMLE<-function(theta)c(F1 = sum((nu_obs-m(theta))/m(theta)),
                            F2 = sum((nu_obs-m(theta))*x/m(theta)))

est_eqMLE2<-function(theta)c(F1 = sum((nu_obs-m2(theta))/m2(theta)),
                            F2 = sum((nu_obs-m2(theta))*(x-break_point)*indicator/m2(theta)))

#Solving the estimating equations
theta_hat<-multiroot(f=est_eqMLE,start=c(5,0),maxiter=1000)$root
theta_hat2<-multiroot(f=est_eqMLE2,start=c(5,0),maxiter=1000)$root

#Fitted mean functions
m_hat_obs<-m(theta_hat)
m2_hat_obs<-m2(theta_hat2)


################################################################################
#Code needed to reproduce Figure 1
################################################################################

par(mar=c(5,5,1,1.5))
plot(x,nu_obs,pch=21,ylab="Photon Counts",xlab="Wavelength",cex.axis=1.5,cex.lab=1.5,type="s",col="grey37",
     ylim=c(0,18))
lines(x,rep(c_hat_obs,K),col="dodgerblue",lwd=3,lty=1)
lines(x,m_hat_obs,col="tomato3",lwd=3,lty=2)
lines(x,m2_hat_obs,col="darkgreen",lwd=3,lty=4)

legend("top",
       legend=c("Constant fit",
                "Linear fit", "Constant+Linear fit"),
       col=c( "dodgerblue", "tomato3","darkgreen"),
       lty=c( 1, 2,4),
       lwd=c( 3, 3,3),
       cex=1.2,horiz = T,
       bty="n") 

################################################################################
#Code needed to reproduce the results reported in Sections 4.1 and 6.3.1
################################################################################

#Constructing Pearson statistic, its linearized version,
#as well as the likelihood ratio (or Cash statistic)
#for testing background uniformity

Pearson_flat <-sum((nu_obs-c_hat_obs)^2/c_hat_obs-1)/sqrt(K)
Pearson_flat

Linearized_flat <-sum((nu_obs-c_hat_obs)/c_hat_obs)/sqrt(K)
Linearized_flat

#Computing the expectation E[nu*log(nu)] at the true value of beta
#and needed to center the likelihood ratio statistic
range<-1:300
Eylny<-sum(range*log(range)*dpois(range, c_hat_obs))
Eylny

Cash_flat<-sum(ifelse(nu_obs>0,nu_obs*log(nu_obs),0)-Eylny-(nu_obs-c_hat_obs)*(log(c_hat_obs)+1))/sqrt(K)

#Simulating the null distribution of the statistics 
#Described above using the parametric bootstrap

#Number of boostrap replicates
B<-100000 #<----100000 in the paper

Pearson0_flat<-Pearson_proj0_flat<-Linearized_flat0<-Cash_flat0<-c()
for(b in 1:B){
  nu0<-rpois(K,c_hat_obs)
  c_hat0<-sum(nu0)/K
  Pearson0_flat[b]<-sum((nu0-c_hat0)^2/c_hat0-1)/sqrt(K)
  Pearson_proj0_flat[b]<-sum((nu0-c_hat_obs)^2/c_hat_obs-1)/sqrt(K)-sum((nu0-c_hat_obs)/c_hat_obs)/sqrt(K)
  Eylny0<-sum(range*log(range)*dpois(range,  c_hat0))
  Cash_flat0[b]<-sum(ifelse(nu0>0,nu0*log(nu0),0)-Eylny0-(nu0-c_hat0)*(log(c_hat0)+1))
  Linearized_flat0[b]<-sum((nu0-c_hat0)/c_hat0)/sqrt(K)
  print(b)
}

########Doing the test using the boostrap and Gaussian approximation for Pearson

pvalue_pearson_flat_boostrap <-2*min(c(mean(Pearson0_flat>Pearson_flat ),mean(Pearson0_flat<Pearson_flat )) )
pvalue_pearson_flat_gaussian <-2*(1-pnorm(abs(Pearson_flat)/sqrt(2)))

pvalue_pearson_flat_boostrap
pvalue_pearson_flat_gaussian

########Doing the test using the boostrap for linearized Pearson statistic and Cash

pvalue_Linearized_flat_boostrap <-2*min(c(mean(Linearized_flat0>Linearized_flat ),mean(Linearized_flat0<Linearized_flat )) )
pvalue_Cash_flat_boostrap <-2*min(c(mean(Cash_flat0>Cash_flat ),mean(Cash_flat0<Cash_flat )) )

pvalue_Linearized_flat_boostrap
pvalue_Cash_flat_boostrap

################################################################################
#Code needed to reproduce the analyses in Section 7.1.1
################################################################################

#Constructing different test statistics, including
#the KS statistics to assess the validity of 
#all the three mean functions considered

Pearson <-sum((nu_obs-m_hat_obs)^2/m_hat_obs-1)/sqrt(K)
Pearson_ps <-max(abs(cumsum((nu_obs-m_hat_obs)^2/m_hat_obs-1)/sqrt(K)))
Linear <-sum((nu_obs-m_hat_obs)/m_hat_obs)/sqrt(K)
Linear_ps <-max(abs(cumsum((nu_obs-m_hat_obs)/m_hat_obs)/sqrt(K)))

Pearson_ps_flat<-max(abs(cumsum((nu_obs-c_hat_obs)^2/c_hat_obs-1)/sqrt(K)))
Linear_ps_flat<-max(abs(cumsum((nu_obs-c_hat_obs)/c_hat_obs)/sqrt(K)))

Pearson_break<-sum((nu_obs-m2_hat_obs)^2/m2_hat_obs-1)/sqrt(K)
Pearson_ps_break<-max(abs(cumsum((nu_obs-m2_hat_obs)^2/m2_hat_obs-1)/sqrt(K)))
Linear_break<-sum((nu_obs-m2_hat_obs)/m2_hat_obs)/sqrt(K)
Linear_ps_break<-max(abs(cumsum((nu_obs-m2_hat_obs)/m2_hat_obs)/sqrt(K)))

############Simulating the null distribution of the test statistics 
############Via the parametric bootstrap

#Function generating observed Poisson counts given a vector of means
rPOI<-function(mui){rpois(1, lambda=mui)}
rPOI<-Vectorize(rPOI,"mui")

#Simulating the null distribution of 
#the Kolmogorov statistics using both Pearson and 
#its linearized counterpart for testing a constant mean function

Pearson_ps0_flat<-c()
Linear_ps0_flat<-c()

for(b in 1:B){#b=1
  nu0<-rPOI(rep(c_hat_obs,K))
  c_hat<-sum(nu0)/K
  Pearson_ps0_flat[b]<-max(abs(cumsum((nu0-c_hat)^2/c_hat-1)/sqrt(K)))
  Linear_ps0_flat[b]<-max(abs(cumsum((nu0-c_hat)/c_hat)/sqrt(K)))
  print(b)
}

#Simulating the null distribution of Pearson and 
#its linearized counterpart as well as the coreesponding KS statistics
#for testing a linear mean function

Pearson0<-Pearson_ps0<-c()
Linear0<-Linear_ps0<-c()

for(b in 1:B){#b=1
  nu0<-rPOI(m_hat_obs)
  est_eqMLE0<-function(theta)c(F1 = sum((nu0-m(theta))/m(theta)),
                               F2 = sum((nu0-m(theta))*x/m(theta)))
  theta_hat0<-multiroot(f=est_eqMLE0,start=c(5,0),maxiter=1000)$root
  m_hat<-m(theta_hat0)
  Pearson0[b]<-sum((nu0-m_hat)^2/m_hat-1)/sqrt(K)
  Pearson_ps0[b]<-max(abs(cumsum((nu0-m_hat)^2/m_hat-1)/sqrt(K)))
  Linear0[b]<-sum((nu0-m_hat)/m_hat)/sqrt(K)
  Linear_ps0[b]<-max(abs(cumsum((nu0-m_hat)/m_hat)/sqrt(K)))
  print(b)
}

#Simulating the null distribution of Pearson and 
#its linearized counterpart as well as the coreesponding KS statistics
#for testing a piecewise linear mean function

Pearson_break0<-Pearson_ps0_break<-c()
Linear_break0<-Linear_ps0_break<-c()

for(b in 1:B){#b=1
  nu0<-rPOI(m2_hat_obs)
  est_eqMLE0<-function(theta)c(F1 = sum((nu0-m2(theta))/m2(theta)),
                               F2 = sum((nu0-m2(theta))*(x-break_point)*indicator/m2(theta)))
  theta_hat0<-multiroot(f=est_eqMLE0,start=c(5,0),maxiter=1000)$root
  m2_hat<-m2(theta_hat0)
  Pearson_break0[b]<-sum((nu0-m2_hat)^2/m2_hat-1)/sqrt(K)
  Pearson_ps0_break[b]<-max(abs(cumsum((nu0-m2_hat)^2/m2_hat-1)/sqrt(K)))
  Linear_break0[b]<-sum((nu0-m2_hat)/m2_hat)/sqrt(K)
  Linear_ps0_break[b]<-max(abs(cumsum((nu0-m2_hat)/m2_hat)/sqrt(K)))
  print(b)
}


#Computing all the p-values 

pvalue_pearson_ps_flat<-mean(Pearson_ps0_flat>Pearson_ps_flat)
pvalue_linear_ps_flat<-mean(Linear_ps0_flat>Linear_ps_flat)

pvalue_pearson <-2*min(c(mean(Pearson0>Pearson ),mean(Pearson0<Pearson )) )
pvalue_linear <-2*min(c(mean(Linear0>Linear ),mean(Linear0<Linear )) )
pvalue_pearson_ps <-mean(Pearson_ps0>Pearson_ps )
pvalue_linear_ps <-mean(Linear_ps0>Linear_ps )

pvalue_pearson_break<-2*min(c(mean(Pearson_break0>Pearson_break ),mean(Pearson_break0<=Pearson_break )) )
pvalue_linear_break<-2*min(c(mean(Linear_break0>Linear_break ),mean(Linear_break0<=Linear_break )) )
pvalue_pearson_ps_break<-mean(Pearson_ps0_break>Pearson_ps_break)
pvalue_linear_ps_break<-mean(Linear_ps0_break>Linear_ps_break)

pvalue_pearson_flat_boostrap
pvalue_Linearized_flat_boostrap 
pvalue_pearson_ps_flat
pvalue_linear_ps_flat

pvalue_pearson 
pvalue_linear 
pvalue_pearson_ps 
pvalue_linear_ps 

pvalue_pearson_break
pvalue_linear_break
pvalue_pearson_ps_break
pvalue_linear_ps_break

#When B=100000

#> pvalue_pearson_flat_boostrap
#[1] 0.57146
#> pvalue_Linearized_flat_boostrap 
#[1] 0.14472
#> pvalue_pearson_ps_flat
#[1] 0.63609
#> pvalue_linear_ps_flat
#[1] 0.00826
#> 
#  > pvalue_pearson 
#[1] 0.48462
#> pvalue_linear 
#[1] 0.23886
#> pvalue_pearson_ps 
#[1] 0.42298
#> pvalue_linear_ps 
#[1] 0.04847
#> 
#  > pvalue_pearson_break
#[1] 0.34104
#> pvalue_linear_break
#[1] 0.99626
#> pvalue_pearson_ps_break
#[1] 0.29373
#> pvalue_linear_ps_break
#[1] 0.42962





