data {
  int<lower=0> n;
  vector[n] x;
  real<lower = 0, upper = 1> t;
  real a;
  real b;
  real pmean;
}

parameters {
  real mu;
  real<lower = 0> sigma_sq;
}

model {
 sigma_sq ~ inv_gamma(a/2, b/2);
 mu ~ normal(pmean, sqrt(sigma_sq));
  target += t * normal_lpdf(x | mu, sqrt(sigma_sq));
}

generated quantities{
  real log_lik;
  log_lik = normal_lpdf(x | mu, sqrt(sigma_sq));
}

