
rm(list = ls())
library(rstan)
#library(shinystan)
#library(spNNGP)       


#------------------------------ NNGP data load ---------------------------------#
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
load(file.choose())


#------------------------------ initialization vals & other vars ---------------#
myinits <-    list(list(beta = c(1, 2, 2), sigma1 = 10, tau = 10, phi1 = 100))
#               list(beta = c(2, 1, 3), sigma1 = 2, tau = 1, phi1 = .2),
#               list(beta = c(0.5, 1.5, .8), sigma1 = 3, tau = 0.5, phi1 = .5))

parameters <- c("beta", "sigmasq1", "tausq", "phi1")


#------------------------------ run stan Variational Bayes ---------------------#
stan_model <- stan_model(file = file.choose())

# ADVI
stan_fit_vb <- vb(stan_model, data=data_nngp, output_samples=1000, seed=123,init=myinits, pars=parameters)


#------------------------------ run MCMC with NUTS -----------------------------#

samples <- stan(
  file = "/Volumes/sghosh/NEcorridor/hierarchicalBayes/nngp_stan/nngp_test.stan",
  data = data_nngp,
  init = myinits,
  pars = parameters,
  iter = 1000,
  chains = 1,
  thin = 100,
  cores=3,
  seed = 123
)




