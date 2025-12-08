library(rstan)
set.seed(2025)

# Harmonic mean estimator

stanfit <- stan("simple model.stan",
                data = list(n = 7,
                            sigma = 2,
                            x = c(0.4, 0.6, 1, -0.3, -0.9, 0.77, 0.12)),
                iter = 10000,
                warmup = 5000)

stanfit



log_lik_sample <- as.matrix(stanfit, pars = "log_lik")
likelihoods <- exp(log_lik_sample)

n <- 7
m <- min(likelihoods)

recips <-  1/likelihoods

# Evidence results

harmonicmean <- 1/mean(recips)

#logharmonicmean <- log(n) - m - log(sum(exp(-log_lik_sample - m)))
#exp(logharmonicmean)

# -----

# Monte Carlo estimator

data <- c(0.4, 0.6, 1, -0.3, -0.9, 0.77, 0.12)
nprior <- 10000

mu_prior <- rnorm(nprior, 0, 1)

likelihood <- numeric(nprior)

for (i in 1:nprior){
  logliki <- sum(dnorm(data, mu_prior[i], 2, log = TRUE))
  likelihood[i] <- exp(logliki)
}

# Evidence results

evidenceMC <- mean(likelihood)
# logevMC <- log(evidenceMC)


# ----

# Analytic calculation

n <- 7
sigmasq <- 4

logev <- -(n/2)*log(8*pi) - 0.5*log((n+4)/4) + ((n*mean(data)/4)^2)/((n+4)/2) - (sum(data^2))/8
