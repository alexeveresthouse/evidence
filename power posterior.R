library(rstan)

# Simple model (Normal with known variance)

# Create vector of temperatures; preferably denser closer to prior
t <- c(0, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
# Create vector to store log-likelihood evaluations in
logliks <- rep(NA, length(t))
# Repeat Stan evaluations at different temperatures
for (i in 1:length(t)){
  posti <- stan("power posterior.stan",
                data = list(n = 7,
                            sigma = 2,
                            t = t[i],
                            x = c(0.4, 0.6, 1, -0.3, -0.9, 0.77, 0.12)))
  
  logliks[i] <- mean(extract(posti)$log_lik)
}

# Update the log evidence using a trapezoidal approximation to the relevant integral
logevsimple <- 0
for (i in 1:(length(t) - 1)){
  logevsimple <- logevsimple + (t[i + 1] - t[i])*(0.5)*(logliks[i + 1] + logliks[i])
}

# Final evidence value
print(exp(logevsimple))


# Normal model with unknown mean and variance

# Create vector of temperatures; preferably denser closer to prior
t <- c(0, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
# Create vector to store log-likelihood evaluations in
logliks <- rep(NA, length(t))
# Repeat Stan evaluations at different temperatures
for (i in 1:length(t)){
  posti <- stan("power posterior complex.stan",
                data = list(n = 7,
                            t = t[i],
                            x = c(0.4, 0.6, 1, -0.3, -0.9, 0.77, 0.12),
                            a = 2,
                            b = 2,
                            pmean = 0),
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
print(exp(logevcomplex))




powerpostcheck

e1 <- mean(extract(powerpostcheck)$log_lik)
e01 <- mean(extract(powerpost0.1)$log_lik)
e025 <- mean(extract(powerpost0.25)$log_lik)
e05 <- mean(extract(powerpost0.5)$log_lik)
e075 <- mean(extract(powerpost0.75)$log_lik)
e09 <- mean(extract(powerpost0.9)$log_lik)

t <- c(0.1, 0.25, 0.5, 0.75, 0.9, 1)
e <- c(e01, e025, e05, e075, e09, e1)

logev <- 0



exp(logev)







library(durhamSLR)
diagnostics(as.array(powerpost))

diagnostics(as.array(powerpostcheck))
diagnostics(as.array(stanfit))
