data {
  int<lower = 0> n;
  int<lower = 1> m;
  matrix[n, m] x;
  array[n] int<lower = 0> y;
  real <lower = 0, upper = 1> t;
}

parameters {
  vector[m] beta;
  real<lower = 0> phi;
}

model {
  // Priors
  vector[n] log_eta = x * beta;
  
  beta ~ normal(0, 1.5);
  phi ~ exponential(0.5);

  // Likelihood
  
  for (i in 1:n){
    target += t * neg_binomial_2_log_lpmf(y[i] | log_eta[i], phi);
  }
}

generated quantities {
  real log_lik = 0;
  vector[n] log_eta = x * beta;
  for (i in 1:n){
      log_lik += neg_binomial_2_log_lpmf(y[i] | log_eta[i], phi);
    }
}
