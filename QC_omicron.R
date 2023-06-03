################################################################################################
# Application to Quebec, Canada in "An SIR-based Bayesian Framework for COVID-19 Infection Estimation"
# Authors: Haoyu Wu, David A. Stephens, Erica E.M. Moodie
################################################################################################

# load in necessary libraries
library(tidyverse)
library(rstan)
options(mc.cores = 3)
rstan_options(auto_write = TRUE)
library(scales)
Sys.setlocale("LC_TIME", "English")

# load in data
load("QC_omicron_data.RData") #load in QC deaths and serosurvey data
qc.omicron <- qc.omicron[qc.omicron$date < as.Date("2022-11-30"), ] #select data between January 01, 2022 and November 30, 2022
n_days <- nrow(qc.omicron)

# population size
N <- 8637650 #Quebec population in T1 2022 (https://statistique.quebec.ca/fr/fichier/bilan-demographique-quebec-edition-2022.pdf)


##------------------------------------------------data pre-processing------------------------------------------------

dates <- as.Date(qc.omicron$date)
n_d_days <- n_days #number of days having death data

qc.omicron$cumdeath <- cumsum(qc.omicron$death)

dt <- qc.omicron$death #daily deaths

# time to death distribution - Gamma
cutoff <- 51
ttd.mean <- 20
ttd.sd <- 10
denom <- pgamma(cutoff, shape = ttd.mean^2/ttd.sd^2, rate = ttd.mean/ttd.sd^2)
ttd <- diff(pgamma(0:cutoff, shape = ttd.mean^2/ttd.sd^2, rate = ttd.mean/ttd.sd^2)) / denom
ttd %*% 1:cutoff
ttd <- c(0, ttd)
tau_rev <- rev(ttd)
plot(ttd)

# add days for inference of the incidence prior to January 01, 2022, which contributed to the deaths that happened at the beginning
add_days <- cutoff

n_d_days <- n_days

n_days <- n_days + add_days
days <- 1:n_days

start_date <- 1
end_date <- n_days + add_days


##------------------------------------------------prior setup------------------------------------------------

# Note: the prior for initial state S_0 and I_0 is directly set in the stan code; the truncations are also in the stan code

# IFR prior
ifr_min <- 0
ifr_max <- 0.04
sd_ifr <- 0.01 
mu_ifr <- 0.01

# AR1 sigma_beta prior
lc_sigma_beta <- 0
sc_sigma_beta <- 1

# AR1 phi_beta prior
phi_beta_min <- 0
phi_beta_max <- 1

# AR1 mu_beta prior
mu_mu_beta <- 0.4
sd_mu_beta <- 0.2
mu_beta_min <- 0

# gamma prior
mu_gamma_inv <- 8
sd_gamma_inv <- 3


##------------------------------------------------serosurveys------------------------------------------------

# number of serosurveys
n_sero <- nrow(qc.omicron.serosurvey)

# serosurvey sampling date
T_s <- qc.omicron.serosurvey$sample_date
t_s <- c()
for (is in 1:n_sero) {
  t_s[is] <- which(as.Date(dates) == as.Date(T_s[is])) + add_days + start_date - 1
}

# serosurvey sample size and positive samples
n_s <- qc.omicron.serosurvey$sample_size
y_s <- qc.omicron.serosurvey$num_positive

##------------------------------------------------model fit------------------------------------------------

# input data
qc_dat <- list(n_days = n_days, 
               n_d_days = n_d_days,
               days = days, 
               n_sero = n_sero,
               N = N, 
               d = dt, 
               n_s = n_s,
               y_s = y_s,
               t_s = t_s,
               #tau
               cutoff = cutoff,
               tau_rev = tau_rev,
               #prior
               ifr_min = ifr_min, 
               ifr_max = ifr_max,
               sd_ifr = sd_ifr,
               mu_ifr = mu_ifr,
               sc_sigma_beta = sc_sigma_beta,
               lc_sigma_beta = lc_sigma_beta,
               phi_beta_min = phi_beta_min,
               phi_beta_max = phi_beta_max,
               sd_mu_beta = sd_mu_beta,
               mu_mu_beta = mu_mu_beta,
               sd_gamma_inv = sd_gamma_inv,
               mu_gamma_inv = mu_gamma_inv)

# initial values
init_fn <- function() {
  list(gamma_inv = 8, beta = rep(0.2, n_days),
       ifr = 0.01)
}

# stan fit
fit <- stan(file = 'QC_omicron_seroonly_ar1_final.stan', init = init_fn,
            data = qc_dat, chains = 3, iter = 30000, warmup = 20000,
            seed = 1898)

# posterior samples
fitdraws <- rstan::extract(fit)

# save posterior results
save(fit, file ="QC_omicron_seroonly_ar1_fit.RData")
save(fitdraws, file ="QC_omicron_seroonly_ar1_fitdraws.RData")


##------------------------------------------------plots------------------------------------------------
load("QC_omicron_seroonly_ar1_fitdraws.RData")

add_days <- 51

# posterior summary
cum_incid_post <- fitdraws$y[,,2] + fitdraws$y[,,3]
cum_incid_med <- apply(cum_incid_post, 2, median)[-(1:add_days)]
cum_incid_cil <- apply(cum_incid_post, 2, quantile, prob = 0.025)[-(1:add_days)]
cum_incid_ciu <- apply(cum_incid_post, 2, quantile, prob = 0.975)[-(1:add_days)]

St_post <- fitdraws$y[,,1]
St_med <- apply(St_post, 2, median)[-(1:add_days)]
St_cil <- apply(St_post, 2, quantile, prob = 0.025)[-(1:add_days)]
St_ciu <- apply(St_post, 2, quantile, prob = 0.975)[-(1:add_days)]

It_post <- fitdraws$y[,,2]
It_med <- apply(It_post, 2, median)[-(1:add_days)]
It_cil <- apply(It_post, 2, quantile, prob = 0.025)[-(1:add_days)]
It_ciu <- apply(It_post, 2, quantile, prob = 0.975)[-(1:add_days)]

Rt_post <- fitdraws$y[,,3]
Rt_med <- apply(Rt_post, 2, median)[-(1:add_days)]
Rt_cil <- apply(Rt_post, 2, quantile, prob = 0.025)[-(1:add_days)]
Rt_ciu <- apply(Rt_post, 2, quantile, prob = 0.975)[-(1:add_days)]

nu_post <- fitdraws$nu
nu_med <- apply(nu_post, 2, median)[-(1:add_days)]
nu_cil <- apply(nu_post, 2, quantile, prob = 0.025)[-(1:add_days)]
nu_ciu <- apply(nu_post, 2, quantile, prob = 0.975)[-(1:add_days)]

beta_post <- fitdraws$beta
beta_med <- apply(beta_post, 2, median)[-(1:add_days)]
beta_cil <- apply(beta_post, 2, quantile, prob = 0.025)[-(1:add_days)]
beta_ciu <- apply(beta_post, 2, quantile, prob = 0.975)[-(1:add_days)]

rt_post <- apply(fitdraws$beta, 2, function(x){x/fitdraws$gamma})
rt_med <- apply(rt_post, 2, median)[-(1:add_days)]
rt_cil <- apply(rt_post, 2, quantile, prob = 0.025)[-(1:add_days)]
rt_ciu <- apply(rt_post, 2, quantile, prob = 0.975)[-(1:add_days)]

qc.omicron.cum.incid <- data.frame("date" = c(as.Date("2022-01-01"), dates), cum_incid_med, cum_incid_cil, cum_incid_ciu)
qc.omicron.St <- data.frame("date" = c(as.Date("2022-01-01"), dates), St_med, St_cil, St_ciu)
qc.omicron.It <- data.frame("date" = c(as.Date("2022-01-01"), dates), It_med, It_cil, It_ciu)
qc.omicron.Rt <- data.frame("date" = c(as.Date("2022-01-01"), dates), Rt_med, Rt_cil, Rt_ciu)
qc.omicron.nu <- data.frame("date" = dates, nu_med, nu_cil, nu_ciu)
qc.omicron.rt <- data.frame("date" = dates, rt_med, rt_cil, rt_ciu)

qc.omicron.St$compartment <- rep("St", nrow(qc.omicron.St))
qc.omicron.It$compartment <- rep("It", nrow(qc.omicron.It))
qc.omicron.Rt$compartment <- rep("Rt", nrow(qc.omicron.Rt))

colnames(qc.omicron.St) <- c("date", "med", "cil", "ciu", "compartment")
colnames(qc.omicron.It) <- c("date", "med", "cil", "ciu", "compartment")
colnames(qc.omicron.Rt) <- c("date", "med", "cil", "ciu", "compartment")

qc.omicron.SIR <- rbind(qc.omicron.St, qc.omicron.It, qc.omicron.Rt)
qc.omicron.SIR$compartment <- factor(qc.omicron.SIR$compartment, levels = c("St", "It", "Rt"), ordered = TRUE)

qc.omicron.St <- data.frame("date" = c(as.Date("2022-01-01"), dates), St_med, St_cil, St_ciu)
qc.omicron.It <- data.frame("date" = c(as.Date("2022-01-01"), dates), It_med, It_cil, It_ciu)
qc.omicron.Rt <- data.frame("date" = c(as.Date("2022-01-01"), dates), Rt_med, Rt_cil, Rt_ciu)

qc.omicron.SIR$serosurvey <- NA
for (i in 1:nrow(qc.omicron.SIR)) {
  for (j in 1:nrow(qc.omicron.serosurvey)) {
    if(as.Date(qc.omicron.SIR$date[i]) == as.Date(qc.omicron.serosurvey$sample_date[j])){
      qc.omicron.SIR$serosurvey[i] <- qc.omicron.serosurvey$seroprevalence[j]
    }
  }
}

# plot

ggplot(data = qc.omicron.nu, aes(x = date, y = nu_med/N*10000, ymax = nu_ciu/N*10000, ymin = nu_cil/N*10000)) + 
  geom_line(col = "#fb8072") +
  geom_ribbon(alpha = 0.2, color = NA, fill = "#fb8072") +
  scale_x_date(date_breaks = "2 month", limits = c(as.Date("2022-01-01"), as.Date("2022-11-30"))) +
  scale_y_continuous(labels = comma, limits = c(0, 160)) +
  xlab("") +
  ylab("Daily new infections per 10,000 people") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  guides(group = "none")


ggplot(data = qc.omicron.cum.incid, aes(x = date, y = cum_incid_med/N, ymax = cum_incid_ciu/N, ymin = cum_incid_cil/N)) + 
  geom_line(col = "#fdb462") +
  geom_ribbon(alpha = 0.2, color = NA, fill = "#fdb462") +
  scale_x_date(date_breaks = "2 month", limits = c(as.Date("2022-01-01"), as.Date("2022-11-30"))) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0,1)) + 
  xlab("") +
  ylab("Cumulative incidence") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  guides(group = "none")


ggplot(data = qc.omicron.SIR, aes(x = date, y = med/N, ymax = ciu/N, ymin = cil/N, col = compartment, fill = compartment)) + 
  geom_line() +
  geom_ribbon(alpha = 0.2, color = NA) +
  geom_point(aes(x = date, y = serosurvey), col = "#33a02c", show.legend = FALSE) +
  scale_x_date(date_breaks = "2 month", limits = c(as.Date("2022-01-01"), as.Date("2022-11-30"))) +
  scale_color_manual(values = c("#80b1d3", "#fdb462", "#8dd3c7", ""), labels = c(bquote(S[t]), bquote(I[t]), bquote(R[t])), name = "") + 
  scale_fill_manual(values = c("#80b1d3", "#fdb462", "#8dd3c7", ""), labels = c(bquote(S[t]), bquote(I[t]), bquote(R[t])), name = "") + 
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0,1)) + 
  xlab("") +
  ylab("Proportion") +
  theme_bw() +
  theme(legend.position = c(0.7, 0.9), legend.direction = "horizontal",
        axis.text.x = element_text(angle = 45, hjust=1)) +
  guides(group = "none")


ggplot(data = qc.omicron.rt, aes(x = date, y = rt_med, ymax = rt_ciu, ymin = rt_cil)) + 
  geom_line(col = "#cab2d6") +
  geom_ribbon(alpha = 0.2, color = NA, fill = "#cab2d6") +
  scale_x_date(date_breaks = "2 month", limits = c(as.Date("2022-01-01"), as.Date("2022-11-30"))) +
  #scale_y_continuous(labels = comma) +
  xlab("") +
  ylab("Reproductive number r(t)") +
  ylim(c(0, 9)) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  guides(group = "none")


