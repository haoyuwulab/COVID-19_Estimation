################################################################################################
# Simulations of "An SIR-based Bayesian Framework for COVID-19 Infection Estimation"
# Authors: Haoyu Wu, David A. Stephens, Erica E.M. Moodie
# Part I: Data generation
################################################################################################

# libraries
library(tidyverse) #expand_grid and ggplot function

##------------------------------------------------data generating parameter setup------------------------------------------------

###------------------------------------------------SEIR parameters------------------------------------------------

# (a) Setting 1: average transmission rate = 0.3, incubation rate = 1/5, number of days: 180

mu_beta <- 0.3
ome <- 1/5
n_days <- 180

# # (b) Setting 2: average transmission rate = 0.3, incubation rate = 1/2, number of days: 120
# 
# mu_beta <- 0.3
# ome <- 1/2
# n_days <- 120
# 
# # (c) Setting 3: average transmission rate = 0.5, incubation rate = 1/5, number of days: 90
# 
# mu_beta <- 0.5
# ome <- 1/5
# n_days <- 90
# 
# # (d) Setting 4: average transmission rate = 0.5, incubation rate = 1/2, number of days: 70
# 
# mu_beta <- 0.5
# ome <- 1/2
# n_days <- 70


# population size
N <- 1000000

# Initial state
S0 <- 0.998*N
E0 <- 0.001*N
I0 <- 0.001*N
R0 <- N - S0 - E0 - I0

# transmission rate beta_t - stationary AR(1)
phi_beta <- 0.9
phi0_beta <- mu_beta * (1 - phi_beta)
sig_beta <- 0.03

be0 <- 1.5 #initial value
be <- c()
be[1] <- be0
nburn_be <- 100 #number of burn-ins
set.seed(18980305) #set seed for reproducibility
for (t in 2:(n_days+nburn_be)) {
  be[t] <- abs(rnorm(1, phi0_beta + phi_beta * be[t-1], sd = sig_beta))
} #randomly generate time series from AR(1)
be <- be[-(1:nburn_be)] #remove burn-in values
plot.ts(be)

# removal rate gamma
gam <- 1/7

###------------------------------------------------death parameters------------------------------------------------

# IFR
IFR <- 0.01

# time to death distribution tau - Gamma distribution with a mean of 15 and a standard deviation of 6
cutoff <- 35
ttd.mean <- 15
ttd.sd <- 6
denom <- pgamma(cutoff, shape = ttd.mean^2/ttd.sd^2, rate = ttd.mean/ttd.sd^2)
ttd <- diff(pgamma((-1):cutoff, shape = ttd.mean^2/ttd.sd^2, rate = ttd.mean/ttd.sd^2)) / denom
ttd %*% 0:cutoff
tau_rev <- rev(ttd)

plot(ttd)

###------------------------------------------------serosurvey parameters------------------------------------------------

# population seroprevalence at the time of the serosurveys: R_{t_s}

# (a) Setting 1: two serosurveys

n_sero <- 2
thetas <- c(0.1, 0.5)

# # (b) Setting 2: five serosurveys
# 
# n_sero <- 5
# thetas <- seq(0.1, 0.7, length.out = 5)

# survey sample size n_s
ns <- c(200, 1000)

# survey bias lambda
pis <- c(0.75, 1, 1.25)

#serosurvey parameter grid
pars <- as.matrix(expand_grid(ns, pis))

##------------------------------------------------transformed parameters------------------------------------------------

# SEIR compartments
St <- Et <- It <- Rt <- rep(0, n_days+1)

# daily incidence nu_t
nut <- rep(0, n_days)

# Poisson convolution mean
lamt <- rep(0, n_days)

# initial state
St[1] <- S0
Et[1] <- E0
It[1] <- I0
Rt[1] <- R0

for (t in 1:n_days) {
  St[t+1] <- St[t] - be[t]*St[t]*It[t]/N
  Et[t+1] <- Et[t] + be[t]*St[t]*It[t]/N - ome*Et[t]
  It[t+1] <- It[t] + ome*Et[t] - gam*It[t]
  Rt[t+1] <- Rt[t] + gam*It[t]
  nut[t] <- ome*Et[t]
  if (t <= cutoff) {
    lamt[t] <- IFR * (t(nut[1:t]) %*% tau_rev[(cutoff-t+2):(cutoff+1)] + It[1] * tau_rev[cutoff-t+1])
  }else{
    lamt[t] <- IFR * (t(nut[(t - cutoff):t]) %*% tau_rev)
  }
}

##------------------------------------------------data pre-processing------------------------------------------------

days <- 1:n_days

start_date <- 1
end_date <- days[length(days)]

##------------------------------------------------prior setup------------------------------------------------

# Note: the prior for initial state S_0 is directly set in the stan code; the truncations are also in the stan code

# (a) Alternative 1: non- or weakly-informative priors

# IFR
ifr_min <- 0
ifr_max <- 0.04

#gamma
gamma_min <- 0
gamma_max <- 1

#AR(1) mu_beta
mu_mu_beta <- 0.5
sd_mu_beta <- 1

# # (b) Alternative 2: informative priors
#  
# # IFR
# ifr_min <- 0
# ifr_max <- 0.04
# sd_ifr <- 0.002 
# mu_ifr <- 0.01
# 
# #gamma
# mu_gamma_inv <- 7
# sd_gamma_inv <- 1.5
# 
# #AR(1) mu_beta
# mu_mu_beta <- 0.3 #change to 0.5 if the data generating parameter is 0.5
# sd_mu_beta <- 0.15
# 

##------------------------------------------------random data generation------------------------------------------------
# number of replications
nrep <- 1000

# number of serosurvey parameter combinations
npar <- nrow(pars)

y_s <- vector("list", nrep)

dt <- vector("list", nrep)

# serosurvey time t_s
t_s <- c()
for (i_sero in 1:n_sero) {
  t_s[i_sero] <- which(Rt >= thetas[i_sero]*N)[1]
}

# data list
dat_list <- vector("list", nrep)

set.seed(2022) #set seed for reproducibility
for (ipar in 1:npar) {
  # assign serosurvey parameters
  ns_temp <- pars[ipar, 1]
  pis_temp <- pars[ipar, 2]
  n_s <- rep(ns_temp, n_sero)
  pi_s <- pis_temp * Rt[t_s] / N
  
  # random data generation
  for (irep in 1:nrep) {
    # number of positive samples in the serosurveys
    y_s[[irep]] <- rbinom(length(t_s), size = n_s, prob = pi_s)
    # daily deaths
    dt[[irep]] <- rpois(n_days, lamt)
    # combine data
    dat_list[[irep]] <- list(y_s = y_s[[irep]], 
                             dt = dt[[irep]])
  }
  
  # save data - alternative 1: non- or weakly-informative priors
  save(be, gam, IFR, 
       mu_beta, phi_beta, sig_beta,
       n_days, days, n_sero, N, 
       St, It, Rt, nut, lamt,
       dat_list,
       n_s, t_s, pi_s,  
       cutoff, tau_rev, 
       ifr_min, ifr_max,
       gamma_min, gamma_max,
       mu_mu_beta, sd_mu_beta, 
       file = paste0("sim_seir_seroonly_ar1_uninfprior_", ipar, "_data_paral.RData"))

  # # save data - alternative 2: informative priors
  # save(be, gam, IFR, 
  #      mu_beta, phi_beta, sig_beta,
  #      n_days, days, n_sero, N, 
  #      St, It, Rt, nut, lamt,
  #      dat_list,
  #      n_s, t_s, pi_s,  
  #      cutoff, tau_rev, 
  #      ifr_min, ifr_max, mu_ifr, sd_ifr,
  #      mu_gamma_inv, sd_gamma_inv,
  #      phi_beta_min, phi_beta_max,
  #      mu_mu_beta, sd_mu_beta, 
  #      lc_sigma_beta, sc_sigma_beta,
  #      file = paste0("sim_seir_seroonly_ar1_", ipar, "_data_paral.RData"))
}  

# plot the SEIR data generating parameters
sim_seir_dgm <- data.frame(St, Et, It, Rt, "nut" = c(0, nut), "date" = 0:n_days)

ggplot(data = sim_seir_dgm, aes(x = date)) + 
  geom_line(aes(y = St/N, color = "1"), size = 0.9) +
  geom_line(aes(y = Et/N, color = "2"), size = 0.9) +
  geom_line(aes(y = It/N, color = "3"), size = 0.9) + 
  geom_line(aes(y = Rt/N, color = "4"), size = 0.9) + 
  scale_color_manual(name = "", values = c("#80b1d3", "#fdb462", "#fb8072", "#8dd3c7"), 
                     labels = c(expression(S[t]), 
                                expression(E[t]), 
                                expression(I[t]), 
                                expression(R[t]))) + 
  xlab("Day") +
  ylab("Proportion") +
  theme_bw() +
  theme(legend.position = "bottom")
