# Harmonic mean estimator

set.seed(2026)

library(rstan)

newcomplexfit <- stan("new complex model.stan",
                      data = list(n = 7,
                                  x = c(0.4, 0.6, 1, -0.3, -0.9, 0.77, 0.12),
                                  pmean = 0,
                                  a = 1,
                                  b = 1),
                      iter = 10000,
                      warmup = 5000)

log_lik_sample_complex <- as.matrix(newcomplexfit, pars = "log_lik")
likelihoods_complex <- exp(log_lik_sample_complex)

recips_complex <-  1/likelihoods_complex

harmonicmean_complex <- 1/mean(recips_complex)
logharmonicmeanans <- log(harmonicmean_complex)

newcomplexfit

library(durhamSLR)

diagnostics(as.array(newcomplexfit))

# Monte Carlo estimator

library(invgamma)

set.seed(2026)

data <- c(0.4, 0.6, 1, -0.3, -0.9, 0.77, 0.12)
n <- length(data)
nprior <- 10000
a <- 1
b <- 1
pmean <- 0

sigma_prior_sq <- rinvgamma(nprior, a, b)
mu_prior <- rnorm(nprior, pmean, sqrt(sigma_prior_sq))
  
likelihood <- numeric(nprior)
logevMC <- numeric(nprior)

for (i in 1:nprior){
  logliki <- sum(dnorm(data, mu_prior[i], sqrt(sigma_prior_sq[i]), log = TRUE))
  likelihood[i] <- exp(logliki)
  logevMC[i] <- logliki
}

#logevMCresult <- mean(logevMC)
evidenceMC <- mean(likelihood)
logevMCans <- log(evidenceMC)

# -----

# Analytic result

biglog <- sum(data^2)/2 - (n*mean(data))^2/(2*n + 2) + b
logevcomplextrue <- -0.5*n*log(2*pi) - 0.5*log(n + 1) + a*log(b) - log(gamma(a)) + log(gamma(0.5*n + a)) - (0.5*n + a)*log(biglog)
exp(logevcomplex)
