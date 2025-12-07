# Harmonic mean estimator

library(rstan)

newcomplexfit <- stan("new complex model.stan",
                      data = list(n = 7,
                                  x = c(0.4, 0.6, 1, -0.3, -0.9, 0.77, 0.12),
                                  pmean = 0,
                                  a = 2,
                                  b = 2),
                      iter = 10000,
                      warmup = 5000)

log_lik_sample_complex <- as.matrix(newcomplexfit, pars = "log_lik")
likelihoods_complex <- exp(log_lik_sample_complex)

recips_complex <-  1/likelihoods_complex

harmonicmean_complex <- 1/mean(recips_complex)

newcomplexfit

library(durhamSLR)

diagnostics(as.array(newcomplexfit))

# Monte Carlo estimator

data <- c(0.4, 0.6, 1, -0.3, -0.9, 0.77, 0.12)
nprior <- 10000
a <- 2
b <- 2
pmean <- 0

sigma_prior_sq <- rinvgamma(nprior, a/2, b/2)
mu_prior <- rnorm(nprior, pmean, 1/sqrt(sigma_prior_sq))
  
likelihood <- numeric(nprior)
logevMC <- numeric(nprior)

for (i in 1:nprior){
  logliki <- sum(dnorm(data, mu_prior[i], sqrt(sigma_prior_sq[i]), log = TRUE))
  likelihood[i] <- exp(logliki)
  logevMC[i] <- logliki
}

logevMCresult <- mean(logevMC)
evidenceMC <- mean(likelihood)
