#This file generates the data for the simulations in Section 4 and the supplement 
#of the paper "Admissible online closed testing must employ e-values"

library(tibble)
library(reshape)
library(matrixStats)

###load functions
source("boosted_e_vals.R")
source("SeqE-Guard.R")

###Set seed to make the results reproducible
set.seed(12345)

m = 1000                #Number of Trials
n = 1000                #Number of Hypotheses per Trial
mu_N=0                  #Mean of the null
mu_As=c(2,3,4)          #Means of the alternative
pi_As=c(0.1,0.3,0.5)    #Proportions of false hypotheses

alpha = 0.1             #Overall significance level

#Parameters for the online simple method
ind_alpha=alpha
a=3
c=log(1/alpha)/(log(1+log(1/alpha)/a)*a)
theta_c=log(1/alpha)/(c*a)

#Parameters for the u-online-simple and m-online-simple method
a_os=(1:n)
weights_os=(6/(pi^2*a_os^2))
alpha_os=alpha*weights_os
c_os=log(1/alpha_os)/(log(1+log(1/alpha_os)/a_os)*a_os)
theta_os=log(1/alpha_os)/(c_os*a_os)

#Parameters for the u-online-Freedman and m-online-Freedman method
a_freed=2^((0:(n-1))/2)
weights_freed=(6/((pi^2+6)*pmax((0:(n-1)),1)^2))
alpha_freed=alpha*weights_freed
kappa_freed=sqrt(2*a_freed*log(1/alpha_freed))+log(1/alpha_freed)/2
lambda_freed=log(1+kappa_freed/a_freed)

#Procedures to compare
lab=c("online-simple", "closed online-simple", "admissible online-simple", "GRO", "hedged GRO", "boosted GRO", "calibrated", "true proportion",
      "m-online-simple", "u-online-simple", "m-online-Freedman", "u-online-Freedman") 

#Final array for true discovery proportions
TDP_array=array(NA, c(n, length(mu_As), length(pi_As), length(lab)), list(1:n,mu_As, pi_As, lab))

for(k in (1:length(mu_As))){
for(l in (1:length(pi_As))){

mu_A=mu_As[k]
pi_A=pi_As[l]

###Generate p-values and e-values
  p = matrix(, nrow = n, ncol = m)
  rejects_os=matrix(, nrow = n, ncol = m)
  e_os = matrix(, nrow = n, ncol = m)
  e_os_admiss = matrix(, nrow = n, ncol = m)
  e_calib = matrix(, nrow = n, ncol = m)
  e_gro = matrix(, nrow = n, ncol = m)
  e_gro_hedged = matrix(, nrow = n, ncol = m)
  
  e_os_average = array(, dim=c(n, n, m))
  e_freed_average = array(, dim=c(n, n, m))
  
  hypo = matrix(, nrow = n, ncol = m)
  for(j in 1:m){
    hypo[, j]=rbinom(n, 1, pi_A)
    X = rnorm(n)
    Z = mu_N * (hypo[, j] - 1) * (-1) + mu_A * hypo[, j] + X
    p[, j] = pnorm(-Z)
    e_gro[,j]=dnorm(Z,mu_A,1)/dnorm(Z,0,1)
  }
  
  #Initialize bounds for different procedures
  bound_os=matrix(, nrow = n, ncol = m )
  bound_e_os=matrix(, nrow = n, ncol = m )
  bound_e_os_admiss=matrix(, nrow = n, ncol = m )
  bound_e_calib=matrix(, nrow = n, ncol = m )
  bound_gro=matrix(, nrow = n, ncol = m )
  bound_gro_hedged=matrix(, nrow = n, ncol = m )
  bound_gro_hedged_boosted=matrix(, nrow = n, ncol = m )
  
  bound_m_os=matrix(, nrow = n, ncol = m )
  bound_u_os=matrix(, nrow = n, ncol = m )
  bound_m_freed=matrix(, nrow = n, ncol = m )
  bound_u_freed=matrix(, nrow = n, ncol = m )
  
  
  for(j in 1:m){
    #Calculate rejection set
    rejects_os[,j]=p[,j]<=ind_alpha
    idx_rejects_os=which(rejects_os[,j])
    
    #Online-simple bound
    for(i in 1:n){
      bound_os[i,j]=sum(rejects_os[1:i,j])-floor(c*(a+i*ind_alpha))
    }
    
    #Closed online-simple bound
    e_os[,j]=exp(theta_c*(rejects_os[,j]-c*ind_alpha))
    bound_e_os[,j]=SeqE_Guard(e_os[,j],idx_rejects_os)
    
    #Admissible online-simple bound
    e_os_admiss[,j]=e_os[,j]/(exp(theta_c*(1-c*ind_alpha))*ind_alpha+exp(-theta_c*c*ind_alpha)*(1-ind_alpha))
    bound_e_os_admiss[,j]=SeqE_Guard(e_os_admiss[,j],idx_rejects_os)
    
    #m-online-simple bound
    for(i in 1:n){
      e_os_average[i, , j]=exp(theta_os[i]*(rejects_os[, j]-c_os[i]*ind_alpha))/(exp(theta_os[i]*(1-c_os[i]*ind_alpha))*ind_alpha+exp(-theta_os[i]*c_os[i]*ind_alpha)*(1-ind_alpha))
    }
    bound_m_os[,j]=SeqE_Guard_average(e_os_average[, , j], weights_os, idx_rejects_os)
    
    #u-online-simple bound
    for(i in 1:n){
      bound_u_os[i,j]=sum(rejects_os[1:i,j])-min(floor(c_os*(a_os+i*ind_alpha)))
    }
    
    #m-online-freedman bound
    for(i in 1:n){
      e_freed_average[i, , j]=exp(lambda_freed[i]*(rejects_os[, j]-ind_alpha))/(exp(lambda_freed[i]*(1-ind_alpha))*ind_alpha+exp(-lambda_freed[i]*ind_alpha)*(1-ind_alpha))
    }
    bound_m_freed[,j]=SeqE_Guard_average(e_freed_average[, , j], weights_freed, idx_rejects_os)
    
    #u-online-freedman bound
    for(i in 1:n){
      epsilon=log((1+pi^2/6)/alpha)+2*log(1+log(max(i*ind_alpha,1),2))
      bound_u_freed[i,j]=sum(rejects_os[1:i,j])-i*ind_alpha-2*sqrt(max(i*ind_alpha,1))*sqrt(epsilon)-epsilon/2
    }
    
    #Bound for calibrated e-values
    x=0.1
    e_calib[,j]=exp(x*qnorm(1-p[,j])-1/2*x^2 )
    bound_e_calib[,j]=boosted_SeqE_Guard(e_calib[,j],idx_rejects_os, x)
    
    #Bound for raw GRO e-values
    bound_gro[,j]=SeqE_Guard(e_gro[,j],idx_rejects_os)
    
    #Bound for hedged GRO e-values
    taus=c(1/2,(cumsum((e_gro[1:(n-1), j]>1))+1/2)/(2:n))
    e_gro_hedged[, j]=taus*e_gro[, j]+(1-taus)
    bound_gro_hedged[,j]=SeqE_Guard(e_gro_hedged[, j],idx_rejects_os)
    
    #Bound for boosted and hedged GRO e-values
    bound_gro_hedged_boosted[,j]=hedged_boosted_SeqE_Guard(e_gro[,j],idx_rejects_os, mu_A)
      
  }
  
#Calculation of true discoveries and rejections
true_discoveries=apply((hypo==1 & rejects_os==1), 2, cumsum)
false_hypotheses=apply(hypo, 2, cumsum)
rejections=apply(rejects_os, 2, cumsum)

#Saving mean TDP bounds in array
TDP_array[1:n,k, l, "online-simple"]=rowMeans(bound_os/pmax(rejections,1))
TDP_array[1:n,k, l, "closed online-simple"]=rowMeans(bound_e_os/pmax(rejections,1))
TDP_array[1:n,k, l, "admissible online-simple"]=rowMeans(bound_e_os_admiss/pmax(rejections,1))
TDP_array[1:n,k, l, "calibrated"]=rowMeans(bound_e_calib/pmax(rejections,1))
TDP_array[1:n,k, l, "GRO"]=rowMeans(bound_gro/pmax(rejections,1))
TDP_array[1:n,k, l, "hedged GRO"]=rowMeans(bound_gro_hedged/pmax(rejections,1))
TDP_array[1:n,k, l, "boosted GRO"]=rowMeans(bound_gro_hedged_boosted/pmax(rejections,1))
TDP_array[1:n,k, l, "m-online-simple"]=rowMeans(bound_m_os/pmax(rejections,1))
TDP_array[1:n,k, l, "u-online-simple"]=rowMeans(bound_u_os/pmax(rejections,1))
TDP_array[1:n,k, l, "m-online-Freedman"]=rowMeans(bound_m_os/pmax(rejections,1))
TDP_array[1:n,k, l, "u-online-Freedman"]=rowMeans(bound_u_os/pmax(rejections,1))
TDP_array[1:n,k, l, "true proportion"]=rowMeans(true_discoveries/pmax(rejections,1))

}
}

#Transform arrays to data frame
TDP_df=as_tibble(melt(TDP_array))
names(TDP_df)=c("index", "mu_A", "pi_A", "Procedure", "TDP")
TDP_df$mu_A=factor(TDP_df$mu_A, levels=mu_As, labels=c("Weak signal", "Medium signal", "Strong signal"))
TDP_df$pi_A=factor(TDP_df$pi_A, levels=pi_As, labels=c("10% non-nulls", "30% non-nulls", "50% non-nulls"))


#Save data frame
TDP_df$TDP[TDP_df$TDP<=0]=0
save(TDP_df, file="results/TD_data.rda")




