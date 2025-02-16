#This file generates the data for the simulations 
#in Section 3.3 of the paper "Admissible online closed testing must employ e-values"

library(tibble)
library(reshape)

###Set seed to make the results reproducible
set.seed(12345)

m = 1000                #Number of Trials
n = 100                 #Number of Hypotheses per Trial
mu_N=0                  #Mean of the null
mu_As=c(2,3,4)          #Means of the alternative
pi_As=c(0.2,0.5,0.8)    #Proportions of false hypotheses

#Final array for estimate of proportion of false hypotheses
taus_array=array(NA, c(n, length(mu_As), length(pi_As)), list(1:n,mu_As, pi_As))

for(k in (1:length(mu_As))){
for(l in (1:length(pi_As))){

mu_A=mu_As[k]
pi_A=pi_As[l]

###Generate p-values and e-values
  taus = matrix(, nrow = n, ncol = m)
  p = matrix(, nrow = n, ncol = m)
  e_gro = matrix(, nrow = n, ncol = m)
  hypo = matrix(, nrow = n, ncol = m)
  for(j in 1:m){
    hypo[, j]=rbinom(n, 1, pi_A)
    X = rnorm(n)
    Z = mu_N * (hypo[, j] - 1) * (-1) + mu_A * hypo[, j] + X
    p[, j] = pnorm(-Z)
    e_gro[,j]=dnorm(Z,mu_A,1)/dnorm(Z,0,1)
    taus[, j]=c(1/2,(cumsum((e_gro[1:(n-1), j]>1))+1/2)/(2:n))
  }
  
taus_array[1:n,k, l]=rowMeans(taus)
}
}

taus_df=as_tibble(melt(taus_array))
names(taus_df)=c("index", "mu_A", "pi_A", "taus")
taus_df$mu_A=factor(taus_df$mu_A, levels=mu_As, labels=c("Weak signal", "Medium signal", "Strong signal"))
taus_df$pi_A_name=factor(taus_df$pi_A, levels=pi_As, labels=c("20% non-nulls", "50% non-nulls", "80% non-nulls"))


#Save data frame
save(taus_df, file="results/TD_data_taus.rda")




