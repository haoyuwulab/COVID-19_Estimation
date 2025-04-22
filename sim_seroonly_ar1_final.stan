data {
  int<lower=1> n_days; // number of follow-up days
  int<lower=1> days[n_days]; // follow-up days
  real<lower=0> N; // population
  int<lower=0> cutoff; // maximum duration from infection to death
  vector<lower=0>[cutoff+1] tau_rev; // time to death distribution
  int<lower=0> d[n_days]; // deaths
  int<lower=1> n_sero; // number of serosurveys
  real<lower=0> ifr_min; 
  real<lower=0> ifr_max;
  real<lower=0> mu_ifr;
  real<lower=0> sd_ifr;
  real<lower=0> mu_gamma_inv;
  real<lower=0> sd_gamma_inv;
  real<lower=0> mu_mu_beta;
  real<lower=0> sd_mu_beta;

  int<lower=1> n_s[n_sero]; // sample sizes of serosurveys
  int<lower=0> y_s[n_sero]; // number of positive samples in the serosurveys
  int<lower=1> t_s[n_sero]; // days at which serosurveys were conducted
}
parameters {
  real<lower=1> gamma_inv;
  vector<lower=0>[n_days] beta;
  real<lower=0> sigma_beta;
  real<lower=0,upper=1> phi_beta;
  real<lower=0> mu_beta;
  real<lower=ifr_min,upper=ifr_max> ifr;
  real y0_tilde;
}
transformed parameters{
  matrix[n_days+1, 3] y; // matrix of SIR compartments
  real<lower=0> phi0_beta = mu_beta * (1 - phi_beta);
  vector<lower=0>[n_days] nu;
  vector<lower=0,upper=1>[n_sero] sp; // population seroprevalence R_t at the time of the serosurveys
  real<lower=0,upper=1> gamma = 1/gamma_inv;
  real<lower=0> death_par[n_days];
  {
    y[1,1] = N*0.99 + 0.01*N*Phi(y0_tilde); // S_0; Phi is the standard normal cumulative distribution function - cdf of y0_tilde follows uniform (0, 1)
    y[1,2] = N - y[1,1]; // I_0
    y[1,3] = 0; // R_0
    for(t in 1:n_days){
      y[t+1,1] = y[t,1] * (1 - beta[t] * y[t,2] / N);
      y[t+1,2] = y[t,2] * (1 + beta[t] * y[t,1] / N - gamma);
      y[t+1,3] = y[t,3] + gamma * y[t,2];
      nu[t] = beta[t] * y[t,1] * y[t,2] / N;
      if(t <= cutoff){
      death_par[t]=ifr*sum(tau_rev[(cutoff-t+2):(cutoff+1)] .*nu[1:t])+ifr * tau_rev[(cutoff-t+1)]*y[1,2];
      }else{
      death_par[t]=ifr*sum(tau_rev .*nu[(t-cutoff):t]);
      }
    }
  }
  
  sp[1:n_sero] = y[t_s[1:n_sero],3]/N;
}
model {
  //priors
  beta[1] ~ normal(mu_beta, sqrt(sigma_beta^2 / (1 - phi_beta^2)));
  beta[2:n_days] ~ normal(phi0_beta + phi_beta * beta[1:(n_days-1)], sigma_beta);
  sigma_beta ~ cauchy(0, 1);
  phi_beta ~ uniform(0,1);
  mu_beta ~ normal(mu_mu_beta, sd_mu_beta);
  gamma_inv ~ gamma(mu_gamma_inv^2/sd_gamma_inv^2,mu_gamma_inv/sd_gamma_inv^2);
  ifr ~ normal(mu_ifr, sd_ifr);
  y0_tilde ~ normal(0,1);

  //likelihoods
  d ~ poisson(death_par[1:n_days]);
  y_s[1:n_sero] ~ binomial(n_s[1:n_sero], sp[1:n_sero]);
}
