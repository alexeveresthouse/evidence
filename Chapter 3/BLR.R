library(rstan)
library(invgamma)
library(durhamSLR)
library(ggplot2)
library(tidyverse)
library(MASS)
library(mvtnorm)
library(GGally)

##### Functions #####

logSumExp <- function(logx){
  c <- max(logx)
  return(c + log(sum(exp(logx - c))))
}

logMeanExp <- function(logx){
  c <- max(logx)
  return(c + log(mean(exp(logx - c))))
}

##### Data setup #####

penguins <- read.csv("penguins.csv") |> na.omit()
penguins_scaled <- as.data.frame(scale(penguins[, 3:6]))
# head(penguins)
# plot(penguins$flipper_len, penguins$body_mass)

# regularfit <- lm(body_mass ~ bill_len + bill_dep + flipper_len,
#                 data = penguins_scaled)

#summary(regularfit)

#### Initial exploration ####

penguins$species <- as.factor(penguins$species)
levels(penguins$species) <- c("A", "C", "G")
pairsdata <- cbind(penguins_scaled, as.factor(penguins$species))
colnames(pairsdata) <- c("bill_len", "bill_dep", "flipper_len", "body_mass", "species")

ggpairs(pairsdata, aes(colour = species))

#### Analytic value of evidence ####

n <- nrow(penguins_scaled)
x <- cbind(1, as.matrix(penguins_scaled[, c("bill_len", "bill_dep", "flipper_len")]))
y <- as.matrix(penguins_scaled$body_mass)
mu0 <- rep(0, 4)
gamma0 <- diag(1/10, nrow = 4, ncol = 4)
a0 <- 2
b0 <- 1


gamman <- t(x)%*%x + gamma0
mun <- solve(gamman)%*%(gamma0%*%mu0 + t(x)%*%y)
an <- a0 + n/2
bn <- b0 + 0.5*(t(y)%*%y + mu0%*%gamma0%*%mu0 - t(mun)%*%gamman%*%mun)

det0 <- det(gamma0)
detn <- det(gamman)

bn <- as.numeric(bn)

log_analytic <- -(n/2) * log(2 * pi) + 
  0.5 * (log(det0) - log(detn)) + 
  a0 * log(b0) - an * log(bn) + 
  lgamma(an) - lgamma(a0)

##### Stan fit #####

linearfit <- stan("BLR model.stan",
                  data = list(n = nrow(penguins_scaled),
                              x = penguins_scaled[, c("bill_len", "bill_dep", "flipper_len")],
                              y = penguins_scaled$body_mass, 
                              m = 3,
                              a = 2,
                              b = 1,
                              sig_alph_sq = 10,
                              sig_beta_sq = 10),
                  iter = 10000,
                  warmup = 5000)

linearfit

diagnostics(as.array(linearfit))

log_lik_sample_complex <- as.matrix(linearfit, pars = "log_lik")
likelihoods_complex <- exp(log_lik_sample_complex)

recips_complex <-  1/likelihoods_complex

# More stable version of log-harmonic mean
log_h_mean <- (log(length(log_lik_sample_complex)) - matrixStats::logSumExp(-log_lik_sample_complex))


##### Logharmonic - Systematic analysis #####

set.seed(2026)

# Compile model once

logharmonic_compiled <- stan_model("BLR model.stan")

# Set up data

n <- nrow(penguins_scaled)
x <- cbind(1, penguins_scaled[, c("bill_len", "bill_dep", "flipper_len")])
y <- penguins_scaled$body_mass
mu0 <- rep(0, 4)
gamma0 <- diag(1/10, nrow = 4, ncol = 4)
a0 <- 2
b0 <- 1

# Set up other variables for analysis to run

niter <- c(100, 1000, 10000, 100000)
nruns <- 15
fits <- list()

logharmonic <- matrix(nrow = nruns, ncol = length(niter))

# Data list

logharmonic_data <- list(n = n,
                         x = x,
                         y = y,
                         m = 4,
                         a = 2, b = 1,
                         beta_precision = gamma0,
                         beta_mean = 0)

# Sample from model

for (i in 1:nruns){
  set.seed(i)
  for (j in 1:length(niter)){
    fits[[j]] <- sampling(logharmonic_compiled,
                          data = logharmonic_data,
                          iter = 150000,
                          warmup = 150000 - niter[j],
                          chains = 1,
                          refresh = 0,
                          show_messages = FALSE)
    
    log_lik_sample_complex <- as.matrix(fits[[j]], pars = "log_lik")
    logsum_recips <- logSumExp(-log_lik_sample_complex)
    
    logharmonic[i, j] <- log(length(log_lik_sample_complex)) - logsum_recips
    
    cat(sprintf("Completed Run %d, Iteration Level %d\n", i, niter[j]))
  }
}

colnames(logharmonic) <- c("100", "1000", "10000", "100000")

palatinate <- rgb(104, 36, 109, maxColorValue = 255)

df.logharmonic <- as.data.frame(logharmonic)
df.logharmonic <- pivot_longer(df.logharmonic, cols = everything(), names_to = "Number of samples", values_to = "Evidence")

ggplot(df.logharmonic, aes(x = `Number of samples`, y = Evidence, fill = "Evidence")) +
  geom_boxplot() +
  ylab("Log evidence") +
  ylim(c(-251, -235)) +
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
  geom_abline(intercept = log_analytic, slope = 0, col = "red")

###### Power Posterior ######

# Create vector of temperatures; preferably denser closer to prior
ntemps <- 25
c <- 5
t <- numeric()
for (j in 1:ntemps){
  t[j] <- (j/ntemps)^c
}
# Create vector to store log-likelihood evaluations in
logliks <- rep(NA, length(t))
# Repeat Stan evaluations at different temperatures
for (i in 1:length(t)){
  posti <- stan("BLR power posterior.stan",
                data = list(n = nrow(penguins_scaled),
                            m = 4,
                            t = t[i],
                            x = cbind(1, penguins_scaled[, c("bill_len", "bill_dep", "flipper_len")]),
                            y = penguins_scaled$body_mass,
                            a = 2,
                            b = 1,
                            beta_mean = 0,
                            beta_precision = diag(1/10, 4, 4)),
                iter = 10000,
                warmup = 5000)
  
  logliks[i] <- mean(extract(posti)$log_lik)
}

# Update the log evidence using a trapezoidal approximation to the relevant integral
logevcomplex <- 0
for (i in 1:(length(t) - 1)){
  logevcomplex <- logevcomplex + (t[i + 1] - t[i])*(0.5)*(logliks[i + 1] + logliks[i])
}

# Final evidence value
logevcomplex

##### Power posterior - Systematic analysis #####

nruns <- 15
niter <- c(100, 1000, 10000, 100000)
powerposterior <- matrix(nrow = nruns, ncol = length(niter))

powerpost_model <- stan_model("BLR power posterior.stan")


for (i in 1:nruns){
  for (j in 1:length(niter)){
    ntemps <- c(10, 20, 30, 40)
    c <- 5
    t <- numeric()
    t <- ((0:ntemps[j]) / ntemps[j])^c
    logliks <- rep(NA, length(t))
    for (k in 1:length(t)){
      posti <- sampling(powerpost_model,
                        data = list(n = nrow(penguins_scaled),
                                    m = 4,
                                    t = t[k],
                                    x = x,
                                    y = y,
                                    a = 2,
                                    b = 1,
                                    beta_mean = 0,
                                    beta_precision = diag(1/10, 4, 4)),
                        iter = niter[j],
                        warmup = 0.5*niter[j],
                        chains = 1,
                        refresh = 0,
                        show_messages = FALSE)
      
      logliks[k] <- mean(rstan::extract(posti)$log_lik)
    }
    logevcomplex <- 0
    for (l in 1:(length(t) - 1)){
      logevcomplex <- logevcomplex + (t[l + 1] - t[l])*(0.5)*(logliks[l + 1] + logliks[l])
    }
    powerposterior[i, j] <- logevcomplex
    cat(sprintf("Completed Run %d, Iteration Level %d\n", i, niter[j]))
  }
}

colnames(powerposterior) <- c("100", "1000", "10000", "100000")

palatinate <- rgb(104, 36, 109, maxColorValue = 255)

df.powerposterior <- as.data.frame(powerposterior)
df.powerposterior <- pivot_longer(df.powerposterior, cols = everything(), names_to = "Number of samples", values_to = "Evidence")

ggplot(df.powerposterior, aes(x = `Number of samples`, y = Evidence, fill = "Evidence")) +
  geom_boxplot() +
  ylab("Log evidence") +
  ylim(c(-255.5, -250)) +
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
  geom_abline(intercept = log_analytic, slope = 0, col = "red")


##### Monte Carlo - systematic analysis #####

n <- nrow(penguins_scaled)
m <- 4
a <- 2
b <- 1
beta_mean <- rep(0, m)
gamma_prior <- diag(1/10, nrow = 4, ncol = 4)

x <- cbind(1, as.matrix(penguins_scaled[, c("bill_len", "bill_dep", "flipper_len")]))
y <- as.matrix(penguins_scaled$body_mass)

nruns <- 10
niter <- c(100, 1000, 10000, 100000)
montecarlo <- matrix(nrow = nruns, ncol = length(niter))

for (k in 1:nruns){
  set.seed(k)
  for (i in 1:length(niter)){
    sigma_prior_sq <- rinvgamma(niter[i], a, b)
    
    likelihood <- numeric(niter[i])
    logev <- numeric(niter[i])
    
    for (j in 1:niter[i]){
      cov_mat_j <- sigma_prior_sq[j] * solve(gamma_prior)
      beta_prior_j <- rmvnorm(1, beta_mean, cov_mat_j)
      
      predictor <- x %*% as.vector(beta_prior_j)
      
      logliki <- sum(dnorm(y, predictor, sqrt(sigma_prior_sq[j]), log = TRUE))
      logev[j] <- logliki
    }
    montecarlo[k,i] <- logMeanExp(logev)
  }
  cat(sprintf("Completed MC Run %d\n", k))
}

colnames(montecarlo) <- c("100", "1000", "10000", "100000")

palatinate <- rgb(104, 36, 109, maxColorValue = 255)

df.montecarlo <- as.data.frame(montecarlo)
df.montecarlo <- pivot_longer(df.montecarlo, cols = everything(), names_to = "Number of samples", values_to = "Evidence")

ggplot(df.montecarlo, aes(x = `Number of samples`, y = Evidence, fill = "Evidence")) +
  geom_boxplot() +
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
  geom_abline(intercept = log_analytic, slope = 0, col = "red")


##### Chib's Estimator #####

# Setup Data and Priors
n <- nrow(penguins_scaled)
m <- 4
x <- cbind(1, as.matrix(penguins_scaled[, c("bill_len", "bill_dep", "flipper_len")]))
y <- as.matrix(penguins_scaled$body_mass)

a <- 2
b <- 1
beta_mean <- rep(0, m)
gamma_prior <- diag(1/10, m, m)

# Pre-compute constant Gibbs parameters
Lambda_n <- t(x) %*% x + gamma_prior
Lambda_n_inv <- solve(Lambda_n)
m_n <- as.vector(Lambda_n_inv %*% (t(x) %*% y + gamma_prior %*% beta_mean))

a_n <- a + (n / 2) + (m / 2) 

# Run the Gibbs sampler
G <- 10000
beta_samples <- matrix(NA, nrow = G, ncol = m)
sigma2_samples <- numeric(G)

sig2_current <- as.numeric(var(y)) # Initial guess

for(g in 1:G) {
  
  cov_beta <- sig2_current * Lambda_n_inv
  beta_current <- rmvnorm(1, m_n, cov_beta)
  beta_samples[g, ] <- beta_current
  
  beta_diff <- as.vector(beta_current) - beta_mean
  err <- y - (x %*% t(beta_current))
  
  b_n <- b + 0.5 * sum(err^2) + 0.5 * as.numeric(t(beta_diff) %*% gamma_prior %*% beta_diff)
  
  sig2_current <- rinvgamma(1, a_n, b_n)
  sigma2_samples[g] <- sig2_current
}

# High posterior density values for calculation
beta_star <- colMeans(beta_samples)
sig2_star <- mean(sigma2_samples)

log_lik_star <- sum(dnorm(y, x %*% beta_star, sqrt(sig2_star), log = TRUE))

log_prior_sig2 <- dinvgamma(sig2_star, a, b, log = TRUE)
log_prior_beta <- dmvnorm(beta_star, mean = beta_mean, sigma = sig2_star * solve(gamma_prior), log = TRUE)
log_prior_star <- log_prior_sig2 + log_prior_beta

log_post_beta <- dmvnorm(beta_star, mean = m_n, sigma = sig2_star * Lambda_n_inv, log = TRUE)

log_post_sig2_evals <- numeric(G)
for(g in 1:G) {
  beta_g <- beta_samples[g, ]
  beta_diff <- beta_g - beta_mean
  err <- y - (x %*% beta_g)
  
  b_n_g <- b + 0.5 * sum(err^2) + 0.5 * as.numeric(t(beta_diff) %*% gamma_prior %*% beta_diff)
  log_post_sig2_evals[g] <- dinvgamma(sig2_star, a_n, b_n_g, log = TRUE)
}

log_post_sig2 <- logMeanExp(log_post_sig2_evals)
log_post_star <- log_post_beta + log_post_sig2

log_evidence_chib <- log_lik_star + log_prior_star - log_post_star

cat(sprintf("Chib's Estimator:    %.3f\n", log_evidence_chib))

##### Chib's estimator - systematic analysis #####

# Setup Data and Priors
n <- nrow(penguins_scaled)
m <- 4
x <- cbind(1, as.matrix(penguins_scaled[, c("bill_len", "bill_dep", "flipper_len")]))
y <- as.matrix(penguins_scaled$body_mass)

a <- 2
b <- 1
beta_mean <- rep(0, m)
gamma_prior <- diag(1/10, m, m)

# Pre-compute constant Gibbs parameters
Lambda_n <- t(x) %*% x + gamma_prior
Lambda_n_inv <- solve(Lambda_n)
m_n <- as.vector(Lambda_n_inv %*% (t(x) %*% y + gamma_prior %*% beta_mean))

a_n <- a + (n / 2) + (m / 2) 

# Run the Gibbs sampler
nruns <- 15
niter <- c(100, 1000, 10000, 100000)
chib <- matrix(nrow = nruns, ncol = length(niter))

for (i in 1:nruns){
  for (j in 1:length(niter)){
    G <- niter[j]
    
    beta_samples <- matrix(NA, nrow = G, ncol = m)
    sigma2_samples <- numeric(G)
    
    sig2_current <- as.numeric(var(y)) # Initial guess
    
    for(g in 1:G) {
      
      cov_beta <- sig2_current * Lambda_n_inv
      beta_current <- rmvnorm(1, m_n, cov_beta)
      beta_samples[g, ] <- beta_current
      
      beta_diff <- as.vector(beta_current) - beta_mean
      err <- y - (x %*% t(beta_current))
      
      b_n <- b + 0.5 * sum(err^2) + 0.5 * as.numeric(t(beta_diff) %*% gamma_prior %*% beta_diff)
      
      sig2_current <- rinvgamma(1, a_n, b_n)
      sigma2_samples[g] <- sig2_current
    }
    
    # High posterior density values for calculation
    beta_star <- colMeans(beta_samples)
    sig2_star <- mean(sigma2_samples)
    
    log_lik_star <- sum(dnorm(y, x %*% beta_star, sqrt(sig2_star), log = TRUE))
    
    log_prior_sig2 <- dinvgamma(sig2_star, a, b, log = TRUE)
    log_prior_beta <- dmvnorm(beta_star, mean = beta_mean, sigma = sig2_star * solve(gamma_prior), log = TRUE)
    log_prior_star <- log_prior_sig2 + log_prior_beta
    
    log_post_beta <- dmvnorm(beta_star, mean = m_n, sigma = sig2_star * Lambda_n_inv, log = TRUE)
    
    log_post_sig2_evals <- numeric(G)
    for(g in 1:G) {
      beta_g <- beta_samples[g, ]
      beta_diff <- beta_g - beta_mean
      err <- y - (x %*% beta_g)
      
      b_n_g <- b + 0.5 * sum(err^2) + 0.5 * as.numeric(t(beta_diff) %*% gamma_prior %*% beta_diff)
      log_post_sig2_evals[g] <- dinvgamma(sig2_star, a_n, b_n_g, log = TRUE)
    }
    
    log_post_sig2 <- logMeanExp(log_post_sig2_evals)
    log_post_star <- log_post_beta + log_post_sig2
    
    log_evidence_chib <- log_lik_star + log_prior_star - log_post_star
    chib[i, j] <- log_evidence_chib
    cat(sprintf("Completed Run %d, Iteration Level %d\n", i, niter[j]))
  }
}

colnames(chib) <- c("100", "1000", "10000", "100000")

palatinate <- rgb(104, 36, 109, maxColorValue = 255)

df.chib <- as.data.frame(chib)
df.chib <- pivot_longer(df.chib, cols = everything(), names_to = "Number of samples", values_to = "Evidence")

ggplot(df.chib, aes(x = `Number of samples`, y = Evidence, fill = "Evidence")) +
  geom_boxplot() +
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
  geom_abline(intercept = log_analytic, slope = 0, col = "red")


##### Direct comparison #####

niter <- 100000
nruns <- 15
direct <- matrix(nrow = nruns, ncol = 4)

# Set up data

n <- nrow(penguins_scaled)
x <- cbind(1, as.matrix(penguins_scaled[, c("bill_len", "bill_dep", "flipper_len")]))
y <- as.vector(penguins_scaled$body_mass)
mu0 <- rep(0, 4)
gamma0 <- diag(1/10, nrow = 4, ncol = 4)
a0 <- 2
b0 <- 1
a <- 2
b <- 1
m <- 4
beta_mean <- rep(0, m)
gamma_prior <- diag(1/10, nrow = 4, ncol = 4)

## SETUP - Harmonic mean ##

fits <- list()

# Compile model once

logharmonic_compiled <- stan_model("BLR model.stan")

# Data list

logharmonic_data <- list(n = n,
                         x = x,
                         y = y,
                         m = 4,
                         a = 2, b = 1,
                         beta_precision = gamma0,
                         beta_mean = 0)

## SETUP - Chib ##

# Pre-compute constant Gibbs parameters
Lambda_n <- t(x) %*% x + gamma_prior
Lambda_n_inv <- solve(Lambda_n)
m_n <- as.vector(Lambda_n_inv %*% (t(x) %*% y + gamma_prior %*% beta_mean))

a_n <- a + (n / 2) + (m / 2) 

# Run the Gibbs sampler
G <- 10000
beta_samples <- matrix(NA, nrow = G, ncol = m)
sigma2_samples <- numeric(G)

sig2_current <- as.numeric(var(y)) # Initial guess

## SETUP - Power Posterior ##

powerpost_model <- stan_model("BLR power posterior.stan")
c <- 5
t <- seq(0, 1, length = 40)^c

## SETUP - Matrices ##

montecarlodc <- matrix(nrow = nruns, ncol = length(niter))
logharmonicdc <- matrix(nrow = nruns, ncol = length(niter))
powerposteriordc <- matrix(nrow = nruns, ncol = length(niter))
chibdc <- matrix(nrow = nruns, ncol = length(niter))

### MONTE CARLO ###

for (k in 1:nruns){
  set.seed(k)
  for (i in 1:length(niter)){
    sigma_prior_sq <- rinvgamma(niter[i], a, b)
    
    likelihood <- numeric(niter[i])
    logev <- numeric(niter[i])
    
    for (j in 1:niter[i]){
      cov_mat_j <- sigma_prior_sq[j] * solve(gamma_prior)
      beta_prior_j <- rmvnorm(1, beta_mean, cov_mat_j)
      
      predictor <- x %*% as.vector(beta_prior_j)
      
      logliki <- sum(dnorm(y, predictor, sqrt(sigma_prior_sq[j]), log = TRUE))
      logev[j] <- logliki
    }
    montecarlodc[k,i] <- logMeanExp(logev)
  }
  cat(sprintf("Completed MC Run %d\n", k))
}

## HARMONIC MEAN ##

for (i in 1:nruns){
  set.seed(i)
  for (j in 1:length(niter)){
    fits[[j]] <- sampling(logharmonic_compiled,
                          data = logharmonic_data,
                          iter = 150000,
                          warmup = 150000 - niter[j],
                          chains = 1,
                          refresh = 0,
                          show_messages = FALSE)
    
    log_lik_sample_complex <- as.matrix(fits[[j]], pars = "log_lik")
    logsum_recips <- logSumExp(-log_lik_sample_complex)
    
    logharmonicdc[i, j] <- log(length(log_lik_sample_complex)) - logsum_recips
    
    cat(sprintf("Completed Run %d, Iteration Level %d\n", i, niter[j]))
  }
}

## POWER POSTERIOR ##

for (i in 1:nruns){
  for (j in 1:length(niter)){
    logliks <- rep(NA, length(t))
    for (k in 1:length(t)){
      posti <- sampling(powerpost_model,
                        data = list(n = nrow(penguins_scaled),
                                    m = 4,
                                    t = t[k],
                                    x = x,
                                    y = y,
                                    a = 2,
                                    b = 1,
                                    beta_mean = 0,
                                    beta_precision = diag(1/10, 4, 4)),
                        iter = niter[j],
                        warmup = 0.5*niter[j],
                        chains = 1,
                        refresh = 0,
                        show_messages = FALSE)
      
      logliks[k] <- mean(rstan::extract(posti)$log_lik)
    }
    logevcomplex <- 0
    for (l in 1:(length(t) - 1)){
      logevcomplex <- logevcomplex + (t[l + 1] - t[l])*(0.5)*(logliks[l + 1] + logliks[l])
    }
    powerposteriordc[i, j] <- logevcomplex
    cat(sprintf("Completed Run %d, Iteration Level %d\n", i, niter[j]))
  }
}

## CHIB ##

for (i in 1:nruns){
  for (j in 1:length(niter)){
    G <- niter[j]
    
    beta_samples <- matrix(NA, nrow = G, ncol = m)
    sigma2_samples <- numeric(G)
    
    sig2_current <- as.numeric(var(y)) # Initial guess
    
    for(g in 1:G) {
      
      cov_beta <- sig2_current * Lambda_n_inv
      beta_current <- rmvnorm(1, m_n, cov_beta)
      beta_samples[g, ] <- beta_current
      
      beta_diff <- as.vector(beta_current) - beta_mean
      err <- y - (x %*% t(beta_current))
      
      b_n <- b + 0.5 * sum(err^2) + 0.5 * as.numeric(t(beta_diff) %*% gamma_prior %*% beta_diff)
      
      sig2_current <- rinvgamma(1, a_n, b_n)
      sigma2_samples[g] <- sig2_current
    }
    
    # High posterior density values for calculation
    beta_star <- colMeans(beta_samples)
    sig2_star <- mean(sigma2_samples)
    
    log_lik_star <- sum(dnorm(y, x %*% beta_star, sqrt(sig2_star), log = TRUE))
    
    log_prior_sig2 <- dinvgamma(sig2_star, a, b, log = TRUE)
    log_prior_beta <- dmvnorm(beta_star, mean = beta_mean, sigma = sig2_star * solve(gamma_prior), log = TRUE)
    log_prior_star <- log_prior_sig2 + log_prior_beta
    
    log_post_beta <- dmvnorm(beta_star, mean = m_n, sigma = sig2_star * Lambda_n_inv, log = TRUE)
    
    log_post_sig2_evals <- numeric(G)
    for(g in 1:G) {
      beta_g <- beta_samples[g, ]
      beta_diff <- beta_g - beta_mean
      err <- y - (x %*% beta_g)
      
      b_n_g <- b + 0.5 * sum(err^2) + 0.5 * as.numeric(t(beta_diff) %*% gamma_prior %*% beta_diff)
      log_post_sig2_evals[g] <- dinvgamma(sig2_star, a_n, b_n_g, log = TRUE)
    }
    
    log_post_sig2 <- logMeanExp(log_post_sig2_evals)
    log_post_star <- log_post_beta + log_post_sig2
    
    log_evidence_chib <- log_lik_star + log_prior_star - log_post_star
    chibdc[i, j] <- log_evidence_chib
    cat(sprintf("Completed Run %d, Iteration Level %d\n", i, niter[j]))
  }
}

### Direct Comparison Matrix ###

direct <- cbind(as.vector(montecarlodc), as.vector(logharmonicdc),
                as.vector(powerposteriordc), as.vector(chibdc))

colnames(direct) <- c("Monte Carlo", "Harmonic Mean", "Power Posterior", "Chib")

### Boxplot ###

palatinate <- rgb(104, 36, 109, maxColorValue = 255)

df.direct <- as.data.frame(direct)
df.direct <- pivot_longer(df.direct, cols = everything(), names_to = "Method", values_to = "Evidence")

ggplot(df.direct, aes(x = Method, y = Evidence, fill = "Evidence")) +
  geom_boxplot() +
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
  geom_abline(intercept = log_analytic, slope = 0, col = "red")