data{
  int<lower = 0> n;
  int<lower = 0> m;
  matrix[n, m] x; // n x m matrix
  vector[n] y;
  real a;
  real b;
  real beta_mean;
  matrix[m, m] beta_precision;
  real<lower = 0, upper = 1> t;
}

parameters{
  vector[m] beta; // Now a vector of length m (number of covariates)
  real<lower = 0> sig_sq;
}

model{
  // Prior
  sig_sq ~ inv_gamma(a, b);
  beta ~ multi_normal_prec(rep_vector(beta_mean, m), 1/sig_sq * beta_precision);
  
  // Likelihood
  target += t * normal_lpdf(y | x*beta, sqrt(sig_sq));
}

generated quantities{
  real log_lik;
  log_lik = normal_lpdf(y | x*beta, sqrt(sig_sq));
}

