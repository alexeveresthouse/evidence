# Harmonic mean estimator

set.seed(2026)

library(rstan)

compiled_model <- stan_model(file = "new complex fit.stan")

set.seed(2026)

niter <- c(100, 1000, 10000, 100000)

fits <- list()
nruns <- 10
logharmonic <- matrix(nrow = nruns, ncol = length(niter))

for (i in 1:nruns){
  set.seed(i)
  for (j in 1:length(niter)){
    fits[[j]] <- stan("new complex model.stan",
                      data = list(n = 7,
                                  x = c(0.4, 0.6, 1, -0.3, -0.9, 0.77, 0.12),
                                  pmean = 0,
                                  a = 1,
                                  b = 1,
                                  kappa = 0.001),
                      iter = 150000,
                      warmup = 150000 - niter[j])
    
    log_lik_sample_complex <- as.matrix(fits[[j]], pars = "log_lik")
    likelihoods_complex <- exp(log_lik_sample_complex)
    
    recips_complex <-  1/likelihoods_complex
    
    harmonicmean_complex <- 1/mean(recips_complex)
    logharmonic[i, j] <- log(harmonicmean_complex)
  }
}

colnames(logharmonic) <- c("100", "1000", "10000", "100000")

boxplot(logharmonic, ylim = c(-12.5, -7.5), xlab = "Number of samples", ylab = "Evidence estimate")
abline(a = logevcomplextruekappa, b = 0, col = "red")

library(ggplot2)
library(tidyverse)

palatinate <- rgb(104, 36, 109, maxColorValue = 255)

df.logharmonic <- as.data.frame(logharmonic)
df.logharmonic <- pivot_longer(df.logharmonic, cols = everything(), names_to = "Number of samples", values_to = "Evidence")

ggplot(df.logharmonic, aes(x = `Number of samples`, y = Evidence, fill = "Evidence")) +
  geom_boxplot() +
  ylim(c(-12.5, -7.5)) +
  ylab("Log evidence") +
  stat_summary(
    geom = "crossbar", 
    width = 0.75, 
    fatten = 0.5, 
    color = "white", 
    fun.data = function(x){ 
      return(c(y = median(x), ymin = median(x), ymax = median(x))) 
    }
  ) +
  scale_fill_manual(values = palatinate) +
  theme_minimal() +
  theme(legend.position = "none", axis.text = element_text(size = 12), axis.title = element_text(size = 16)) +
  geom_abline(intercept = logevcomplextruekappa, slope = 0, col = "red")

plot(niter, logharmonic, type = "l")

newcomplexfit <- stan("new complex model.stan",
                      data = list(n = 7,
                                  x = c(0.4, 0.6, 1, -0.3, -0.9, 0.77, 0.12),
                                  pmean = 0,
                                  a = 1,
                                  b = 1,
                                  kappa = 0.001),
                      iter = 20000,
                      warmup = 10000)

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
kappa <- 0.001

nruns <- 10

sigma_prior_sq <- rinvgamma(nprior, a, b)
mu_prior <- rnorm(nprior, pmean, sqrt(sigma_prior_sq/kappa))

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

set.seed(2026)

data <- c(0.4, 0.6, 1, -0.3, -0.9, 0.77, 0.12)
n <- length(data)
nprior <- 10000
niter <- c(100, 1000, 10000, 100000, 1000000)
a <- 1
b <- 1
pmean <- 0
kappa <- 0.001

nruns <- 10

montecarlo <- matrix(nrow = nruns, ncol = length(niter))

for (k in 1:nruns){
  set.seed(k)
  for (i in 1:length(niter)){
    sigma_prior_sq <- rinvgamma(niter[i], a, b)
    mu_prior <- rnorm(niter[i], pmean, sqrt(sigma_prior_sq/kappa))
    
    likelihood <- numeric(niter[i])
    logevMC <- numeric(niter[i])
    
    for (j in 1:length(niter)){
      logliki <- sum(dnorm(data, mu_prior[j], sqrt(sigma_prior_sq[j]), log = TRUE))
      likelihood[j] <- exp(logliki)
      logevMC[j] <- logliki
    }
    evidenceMC <- mean(likelihood)
    montecarlo[i,j] <- log(evidenceMC)
  }
}

colnames(montecarlo) <- c("100", "1000", "10000", "100000")

boxplot(montecarlo, xlab = "Number of samples", ylab = "Evidence estimate")
abline(a = logevcomplextruekappa, b = 0, col = "red")


### Gemini ###

# Assuming data, a, b, pmean, kappa, nruns are already defined

nruns <- 100
niter <- c(100, 1000, 10000, 100000) # Example sample sizes
montecarlo <- matrix(nrow = nruns, ncol = length(niter))

for (k in 1:nruns) {
  set.seed(k)
  
  for (i in 1:length(niter)) {
    # 1. Generate N samples for this specific iteration
    current_n <- niter[i]
    sig_sq <- rinvgamma(current_n, a, b)
    mu <- rnorm(current_n, pmean, sqrt(sig_sq/kappa))
    
    # 2. Vectorized likelihood calculation (much faster than a third loop!)
    # For each pair of (mu, sig_sq), calculate the joint likelihood of the data
    lik_samples <- sapply(1:current_n, function(idx) {
      exp(sum(dnorm(data, mean = mu[idx], sd = sqrt(sig_sq[idx]), log = TRUE)))
    })
    
    # 3. Store the log of the mean likelihood (the log-evidence)
    montecarlo[k, i] <- log(mean(lik_samples))
  }
}




# -----

# Analytic result

biglog <- sum(data^2)/2 - (n*mean(data))^2/(2*n + 2) + b
logevcomplextrue <- -0.5*n*log(2*pi) - 0.5*log(n + 1) + a*log(b) - log(gamma(a)) + log(gamma(0.5*n + a)) - (0.5*n + a)*log(biglog)
exp(logevcomplex)

# Analytic result with kappa term

kappa <- 0.001
biglogkappa <- sum(data^2)/2 - (n*mean(data))^2/(2*n + 2*kappa) + b
logevcomplextruekappa <- -0.5*n*log(2*pi) + 0.5*log(kappa) - 0.5*log(n + kappa) + a*log(b) - log(gamma(a)) + log(gamma(0.5*n + a)) - (0.5*n + a)*log(biglogkappa)

