data {
  int<lower = 0> n;
  int<lower = 1> m;
  matrix[n, m] x;
  array[n] int<lower = 0> y;
  real <lower = 0, upper = 1> t;
}

parameters {
  vector[m] beta;
  real logit_p;
}

model {
  // Priors
  real p = inv_logit(logit_p);
  vector[n] log_eta = x * beta;
  
  beta ~ normal(0, 1.5);
  logit_p ~ logistic(0, 1);

  // Likelihood
  
  for (i in 1:n){
    real log_lik_i;
    
    if (y[i] == 0){
      log_lik_i = log_sum_exp(log(p),
      log1m(p) + poisson_log_lpmf(y[i] | log_eta[i]));
    }
    else {
      log_lik_i = log1m(p) + poisson_log_lpmf(y[i] | log_eta[i]);
    }
    
    target += t * log_lik_i;
  }
}

generated quantities {
  real log_lik = 0;
  real p = inv_logit(logit_p);
  vector[n] log_eta = x * beta;
  
  for (i in 1:n){
    if (y[i] == 0) {
      log_lik += log_sum_exp(log(p),
      log1m(p) + poisson_log_lpmf(y[i] | log_eta[i]));
    }
    else {
      log_lik += log1m(p) + poisson_log_lpmf(y[i] | log_eta[i]);
    }
  }
}

