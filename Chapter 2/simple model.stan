data {
  int<lower=0> n;
  vector[n] x;
  real<lower=0> sigma;
}

parameters {
  real mu;
}

model {
  mu ~ normal(0,1);
  x ~ normal(mu, sigma);
}

generated quantities{
  real log_lik;
  log_lik = normal_lpdf(x | mu, sigma);
}

