rm(list = ls())

library(RRDecon)
library(Rcpp)
library(genlasso)

# Some parameter settings for the mixture model that will generate signals
parlist = list()
parlist[[1]] =  list(weights=c(0.5,0.4,.1), mu=c(-1.5,1.5,4), tau2 = c(1,2,2))
parlist[[2]] = list(weights=c(1/3,1/3,1/3), mu=c(0,-2,3), tau2 = c(2,.1,.4))
parlist[[3]] = list(weights=c(0.3,0.4,0.3), mu=-c(0,0,0), tau2 = c(0.1,1,9))
parlist[[4]] = list(weights=c(0.2,0.3,0.3,0.2), mu=c(-3,-1.5,1.5,3), tau2 = c(0.01,.01,.01,0.01))
parlist[[5]] =   list(weights=c(0.4,0.4), mu=c(-1.5,1.5), tau2 = c(1,1))
parlist[[6]] =  list(weights=c(0.5,0.4,0.1), mu=c(0,-2,3), tau2 = c(.2,.1,.4))

bound = c(20,20,20,5,5)
# Sim settings
n_grid = c(100000, 50000, 25000,10000,2000 )
d_grid = c(150,200,250)
NMC = 50
Num_Methods = 5
Num_densities = 4


wasserstein_dist = array(0,  c(length(d_grid),length(n_grid), Num_densities,Num_Methods,NMC))


MSE_out95 = array(0,  c(length(d_grid),length(n_grid), Num_densities,Num_Methods,NMC))
MSE_out99 = array(0,  c(length(d_grid),length(n_grid), Num_densities,Num_Methods,NMC))

KL_out = array(0,  c(length(d_grid),length(n_grid), Num_densities,Num_Methods,NMC))
bayes_out = array(0,  c(length(d_grid),length(n_grid), Num_densities,Num_Methods +3 ,NMC))

efron_mse_95 = array(0,  c(length(d_grid),length(n_grid), Num_densities,NMC))
efron_mse_99 = array(0,  c(length(d_grid),length(n_grid), Num_densities,NMC))
efron_bayes = array(0,  c(length(d_grid),length(n_grid), Num_densities,NMC))
efron_wasserstein_dist   = array(0,  c(length(d_grid),length(n_grid), Num_densities,NMC))

d_ind = 1

n_ind = 1
den = 2
trial = 1


n = n_grid[n_ind]
d = 150



mypars = parlist[[den]]
true_mu = mypars$mu
true_weights =  mypars$weights
true_tau2 = mypars$tau2

mu = rnormix(n,true_weights,true_mu,true_tau2)
hist(mu)

n = length(mu)
y = mu + rnorm(n)


system.time({temp = L2_deconvolution(y,prop = 0.25)})
loc =  temp$loc
f_true =  dnormix(loc,true_weights,true_mu,true_tau2)
mean((f_true -  temp$f_hat)^2)


temp3 = L2_deconvolution_path(y, d= 200)
matplot(t(temp3$y_hat),type = "l")
matplot(t(temp3$f_hat),type ="l")

system.time({temp4 = L1_deconvolution_path(y, d= 150)})
matplot(t(temp4$y_hat),type = "l")
matplot(t(temp4$f_hat),type = "l")


temp5 = L1_deconvolution(y, d= 150)
loc2 =  temp5$loc
f_true2 =  dnormix(loc2,true_weights,true_mu,true_tau2)
mean((f_true2 -  temp5$f_hat)^2)
plot(temp5$f_hat)
