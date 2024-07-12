#This file generates the data for the simulations in Section 4 of the paper "Online true discovery guarantee with e-values"

library(tibble)
library(reshape)

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
a=1
c=log(1/alpha)/(log(1+log(1/alpha)/a)*a)
theta_c=log(1/alpha)/(c*a)

#Procedures to compare
lab=c("online-simple", "closed online-simple", "admissible online-simple", "GRO", "weighted GRO", "boosted GRO", "calibrated", "true proportion") 

#Initialize matrices for true discovery proportions
prop_os=matrix(, nrow=length(mu_As)*length(pi_As), ncol=n)
prop_e_os=matrix(, nrow=length(mu_As)*length(pi_As), ncol=n)
prop_e_os_admiss=matrix(, nrow=length(mu_As)*length(pi_As), ncol=n)
prop_e_calib=matrix(, nrow=length(mu_As)*length(pi_As), ncol=n)
prop_gro=matrix(, nrow=length(mu_As)*length(pi_As), ncol=n)
prop_gro_weighted=matrix(, nrow=length(mu_As)*length(pi_As), ncol=n)
prop_gro_boosted=matrix(, nrow=length(mu_As)*length(pi_As), ncol=n)
prop_gro_weighted_boosted=matrix(, nrow=length(mu_As)*length(pi_As), ncol=n)

prop_true=matrix(, nrow=length(mu_As)*length(pi_As), ncol=n)

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
  e_gro_weighted = matrix(, nrow = n, ncol = m)
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
  bound_gro_weighted=matrix(, nrow = n, ncol = m )
  bound_gro_weighted_boosted=matrix(, nrow = n, ncol = m )
  bound_gro_boosted=matrix(, nrow = n, ncol = m )
  
  
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
    
    #Bound for calibrated e-values
    x=0.1
    e_calib[,j]=exp(x*qnorm(1-p[,j])-1/2*x^2 )
    bound_e_calib[,j]=boosted_SeqE_Guard(e_calib[,j],idx_rejects_os, x)
    
    #Bound for raw GRO e-values
    bound_gro[,j]=SeqE_Guard(e_gro[,j],idx_rejects_os)
    
    #Bound for weighted GRO e-values
    taus=c(1/2,(cumsum((e_gro[1:(n-1), j]>1))+1/2)/(2:n))
    e_gro_weighted[, j]=taus*e_gro[, j]+(1-taus)
    bound_gro_weighted[,j]=SeqE_Guard(e_gro_weighted[, j],idx_rejects_os)
    
    #Bound for boosted GRO e-values
    bound_gro_boosted[,j]=boosted_SeqE_Guard(e_gro[,j],idx_rejects_os, mu_A)
    
    #Bound for boosted and weighted GRO e-values
    bound_gro_weighted_boosted[,j]=weighted_boosted_SeqE_Guard(e_gro[,j],idx_rejects_os, mu_A)
      
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
TDP_array[1:n,k, l, "weighted GRO"]=rowMeans(bound_gro_weighted/pmax(rejections,1))
TDP_array[1:n,k, l, "boosted GRO"]=rowMeans(bound_gro_weighted_boosted/pmax(rejections,1))
TDP_array[1:n,k, l, "true proportion"]=rowMeans(true_discoveries/pmax(rejections,1))

}
}

#Transform array to data frame
TDP_df=as_tibble(melt(TDP_array))
names(TDP_df)=c("index", "mu_A", "pi_A", "Procedure", "TDP")
TDP_df$mu_A=factor(TDP_df$mu_A, levels=mu_As, labels=c("Weak signal", "Medium signal", "Strong signal"))
TDP_df$pi_A=factor(TDP_df$pi_A, levels=pi_As, labels=c("10% non-nulls", "30% non-nulls", "50% non-nulls"))

#Save data frame
save(TDP_df, file="results/TD_data.rda")




