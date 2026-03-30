library(rstan)
library(ggplot2)
library(tidyr)
set.seed(2026)
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
set.seed(2026)

# Create vector of temperatures; preferably denser closer to prior
ntemps <- 50
c <- 5
t <- numeric()
for (x in 1:ntemps){
  t[x] <- (x/ntemps)^c
}
# Create vector to store log-likelihood evaluations in
logliks <- rep(NA, length(t))
# Repeat Stan evaluations at different temperatures
for (i in 1:length(t)){
  posti <- stan("power posterior complex.stan",
                data = list(n = 7,
                            t = t[i],
                            x = c(0.4, 0.6, 1, -0.3, -0.9, 0.77, 0.12),
                            a = 1,
                            b = 1,
                            pmean = 0,
                            kappa = 0.001),
                iter = 10000,
                warmup = 5000)
  
  logliks[i] <- mean(rstan::extract(posti)$log_lik)
}

# Update the log evidence using a trapezoidal approximation to the relevant integral
logevcomplex <- 0
for (i in 1:(length(t) - 1)){
  logevcomplex <- logevcomplex + (t[i + 1] - t[i])*(0.5)*(logliks[i + 1] + logliks[i])
}

# Final evidence value
logevcomplex
print(exp(logevcomplex))

# Analytic result with kappa term

data <- c(0.4, 0.6, 1, -0.3, -0.9, 0.77, 0.12)
n <- 7
a <- 1
b <- 1
kappa <- 0.001
biglogkappa <- sum(data^2)/2 - (n*mean(data))^2/(2*n + 2*kappa) + b
logevcomplextruekappa <- -0.5*n*log(2*pi) + 0.5*log(kappa) - 0.5*log(n + kappa) + a*log(b) - log(gamma(a)) + log(gamma(0.5*n + a)) - (0.5*n + a)*log(biglogkappa)



#### Systematic analysis ####

nruns <- 100
niter <- c(100, 1000, 10000, 100000)
powerposterior <- matrix(nrow = nruns, ncol = length(niter))

for (i in 1:nruns){
  for (j in 1:length(niter)){
    ntemps <- 50
    c <- 5
    t <- numeric()
    for (x in 1:ntemps){
      t[x] <- (x/ntemps)^c
    }
    logliks <- rep(NA, length(t))
    for (k in 1:length(t)){
      posti <- stan("power posterior complex.stan",
                    data = list(n = 7,
                                t = t[k],
                                x = c(0.4, 0.6, 1, -0.3, -0.9, 0.77, 0.12),
                                a = 1,
                                b = 1,
                                pmean = 0,
                                kappa = 0.001),
                    iter = niter[j],
                    warmup = 0.5*niter[j])
      
      logliks[k] <- mean(extract(posti)$log_lik)
    }
    logevcomplex <- 0
    for (l in 1:(length(t) - 1)){
      logevcomplex <- logevcomplex + (t[l + 1] - t[l])*(0.5)*(logliks[l + 1] + logliks[l])
    }
    powerposterior[i, j] <- logevcomplex
  }
}

colnames(powerposterior) <- c("100", "1000", "10000", "100000")

colnames(powerposterior) <- c("100", "1000", "10000", "100000")

palatinate <- rgb(104, 36, 109, maxColorValue = 255)

df.powerposterior <- as.data.frame(powerposterior)
df.powerposterior <- pivot_longer(df.powerposterior, cols = everything(), names_to = "Number of samples", values_to = "Evidence")

ggplot(df.powerposterior, aes(x = `Number of samples`, y = Evidence, fill = "Evidence")) +
  geom_boxplot() +
  ylim(c(-12.5, -12)) +
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

library(MASS)

write.matrix(powerposterior, "powerposterior matrix")
