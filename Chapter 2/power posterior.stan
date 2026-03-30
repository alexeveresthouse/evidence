data {
  int<lower=0> n;
  vector[n] x;
  real<lower=0> sigma;
  real<lower=0, upper=1> t;
}

parameters {
  real mu;
}

model {
  //Prior
  mu ~ normal(0,1);

  // Likelihood with power posterior
  target+= t * normal_lpdf(x | mu, sigma);
}

generated quantities{
  real log_lik;
  log_lik = normal_lpdf(x | mu, sigma);
}