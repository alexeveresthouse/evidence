data {
  int<lower=0> n;
  vector[n] x;
}

parameters {
  real mu;
  real<lower = 0> sigma;
}

model {
  mu ~ normal(0,1);
  sigma ~ inv_gamma(1,2);
  x ~ normal(mu, sigma);
}

generated quantities{
  real log_lik;
  log_lik = normal_lpdf(x | mu, sigma);
}

