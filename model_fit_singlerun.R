################################################################################################
# Simulations of "An SIR-based Bayesian Framework for COVID-19 Infection Estimation"
# Authors: Haoyu Wu, David A. Stephens, Erica E.M. Moodie
# Part II: Model fitting
################################################################################################

# libraries
library(rstan)

# specify the serosurvey parameter combination
ipar <- 1 #choose from 1,2,3,4,5,6;

## 1 - survey sample size = 200; bias parameter lambda = 0.75
## 2 - survey sample size = 200; bias parameter lambda = 1
## 3 - survey sample size = 200; bias parameter lambda = 1.25
## 4 - survey sample size = 1000; bias parameter lambda = 0.75
## 5 - survey sample size = 1000; bias parameter lambda = 1
## 6 - survey sample size = 1000; bias parameter lambda = 1.25

# alternative 1: non- or weakly-informative priors
load(paste0("sim_seir_seroonly_ar1_uninfprior_", ipar, "_data_paral.RData"))

# # alternative 2: informative priors
# load(paste0("sim_seir_seroonly_ar1_", ipar, "_data_paral.RData"))

# number of replications per scenario
nrep <- 1000

# specify which replicated dataset to use - choose from 1 to 1000
irep <- 1 

# set seed for reproducibility
set.seed(141)

# function to generate the initial value for IFR
init_fn <- function() {
  list(ifr=min(max(rnorm(1, 0.015, 0.005), 0.003), 0.035))
}

##------------------------------------------------model fit------------------------------------------------

# alternative 1: non- or weakly-informative priors
# data input into stan
sim_dat <- list(n_days = n_days, 
                days = days, 
                n_sero = n_sero,
                N = N, 
                d = dat_list[[irep]]$dt, 
                n_s = n_s,
                y_s = dat_list[[irep]]$y_s,
                t_s = t_s,
                # tau
                cutoff = cutoff,
                tau_rev = tau_rev,
                # prior
                ifr_min = ifr_min, 
                ifr_max = ifr_max,
                gamma_min = gamma_min,
                gamma_max = gamma_max,
                sd_mu_beta = sd_mu_beta,
                mu_mu_beta = mu_mu_beta)

# fit the model with stan
fit <- stan(file = 'sim_seroonly_ar1_uninfprior_final.stan', init = init_fn,
            data = sim_dat, chains = 3, iter = 10000, warmup = 5000, cores = 3)

# extract posterior samples and summary statistics
fitdraws <- as.data.frame(fit)
fit.sp <- get_sampler_params(fit, inc_warmup = FALSE)
fitdraws$divergent <- c(sapply(fit.sp, FUN = function(x) x[ ,"divergent__"]))
fit.summary <- rstan::summary(fit)

save(fitdraws, file = paste0("sim_seir_seroonly_ar1_uninfprior_", ipar, "_", irep, "_fitdraws.RData"))
save(fit.summary, file = paste0("sim_seir_seroonly_ar1_uninfprior_", ipar, "_", irep, "_summary.RData"))

# # alternative 2: informative priors
# # data input into stan
# sim_dat <- list(n_days = n_days, 
#                 days = days, 
#                 n_sero = n_sero,
#                 N = N, 
#                 d = dat_list[[irep]]$dt, 
#                 n_s = n_s,
#                 y_s = dat_list[[irep]]$y_s,
#                 t_s = t_s,
#                 # tau
#                 cutoff = cutoff,
#                 tau_rev = tau_rev,
#                 # prior
#                 ifr_min = ifr_min, 
#                 ifr_max = ifr_max,
#                 sd_ifr = sd_ifr,
#                 mu_ifr = mu_ifr,
#                 sd_gamma_inv = sd_gamma_inv,
#                 mu_gamma_inv = mu_gamma_inv,
#                 sd_mu_beta = sd_mu_beta,
#                 mu_mu_beta = mu_mu_beta)
#
# # fit the model with stan
# fit <- stan(file = 'sim_seroonly_ar1_final.stan', init = init_fn,
#             data = sim_dat, chains = 3, iter = 10000, warmup = 5000, cores = 3)
#
# # extract posterior samples and summary
# fitdraws <- as.data.frame(fit)
# fit.sp <- get_sampler_params(fit, inc_warmup = FALSE)
# fitdraws$divergent <- c(sapply(fit.sp, FUN = function(x) x[ ,"divergent__"]))
# fit.summary <- rstan::summary(fit)
# 
# save(fitdraws, file = paste0("sim_seir_seroonly_ar1_", ipar, "_", irep, "_fitdraws.RData"))
# save(fit.summary, file = paste0("sim_seir_seroonly_ar1_", ipar, "_", irep, "_summary.RData"))

